"""3D modeler window: renders the project geometry with PyVista and supports click-to-select."""
from __future__ import annotations

import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from PySide6.QtCore import Signal

from .model import AXES, Project, _in_plane_axes

BACKGROUND_TOP = "#dfe7f2"
BACKGROUND_BOTTOM = "#ffffff"
SELECT_COLOR = "#ff00ff"
PORT_COLOR = "#e02020"


def primitive_mesh(obj, variables, project=None) -> pv.PolyData:
    if obj.kind == "cad":
        from . import cad
        pts, tris = cad.mesh(project, obj)
        if len(tris) == 0:
            raise ValueError("no surface to display")
        return pv.PolyData(pts, np.hstack([np.full((len(tris), 1), 3), tris]).ravel())
    v = obj.values(variables)
    c = np.array([v["x"], v["y"], v["z"]])
    a = AXES.index(obj.axis)
    direction = np.eye(3)[a]
    if obj.kind == "box":
        lo, hi = obj.bounds(variables)
        return pv.Box(bounds=(lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]))
    if obj.kind == "cylinder":
        h = v["height"]
        return pv.Cylinder(center=c + direction * h / 2, direction=direction * np.sign(h or 1),
                           radius=v["radius"], height=abs(h), resolution=48)
    if obj.kind == "sphere":
        return pv.Sphere(radius=v["radius"], center=c, theta_resolution=40, phi_resolution=24)
    # rectangle (sheet)
    u, w = _in_plane_axes(a)
    du, dw = np.zeros(3), np.zeros(3)
    du[u], dw[w] = v["w"], v["l"]
    pts = np.array([c, c + du, c + du + dw, c + dw])
    return pv.PolyData(pts, faces=[4, 0, 1, 2, 3])


def port_mesh(port, variables, size_hint=1.0) -> pv.PolyData:
    a, b = (np.array(p) for p in port.corners(variables))
    ext = np.abs(b - a)
    if np.count_nonzero(ext) <= 1:
        # line port: draw as a tube scaled to the model so it stays visible and pickable
        line = pv.Line(a, b)
        return line.tube(radius=max(0.006 * size_hint, 1e-9), n_sides=16)
    return pv.Box(bounds=(min(a[0], b[0]), max(a[0], b[0]), min(a[1], b[1]),
                          max(a[1], b[1]), min(a[2], b[2]), max(a[2], b[2])))


class ModelerView(QtInteractor):
    """PyVista Qt widget showing the model. Emits `picked(name)` on a click (not a drag)."""

    picked = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.set_background(BACKGROUND_BOTTOM, top=BACKGROUND_TOP)
        self.enable_anti_aliasing("ssaa")
        self.add_axes(interactive=False, line_width=2, color="black")
        self._actors: dict[str, object] = {}
        self._actor_names: dict[int, str] = {}
        self._tessellated: dict[str, bool] = {}
        self._colors: dict[str, tuple] = {}
        self._sheet_actors: list = []
        self._selected: str | None = None
        self._show_grid = True
        self._press_pos = None
        self.iren.add_observer("LeftButtonPressEvent", self._on_press)
        self.iren.add_observer("LeftButtonReleaseEvent", self._on_release)

    # ------------------------------------------------------------------ picking
    def _on_press(self, obj, _event):
        self._press_pos = self.iren.get_event_position()

    def _on_release(self, obj, _event):
        pos = self.iren.get_event_position()
        if self._press_pos is None or abs(pos[0] - self._press_pos[0]) + abs(pos[1] - self._press_pos[1]) > 3:
            return
        actor = self._pick_actor(*pos)
        name = self._actor_names.get(id(actor), "") if actor is not None else ""
        self.picked.emit(name)

    def _pick_actor(self, x, y):
        import vtk
        # Ray-cast rather than vtkPropPicker: its hardware selection finds nothing under macOS Cocoa.
        picker = vtk.vtkCellPicker()
        picker.SetTolerance(0.0005)
        picker.Pick(x, y, 0, self.renderer)
        actor = picker.GetActor()
        if actor is None or not self._sheet_actors:
            return actor
        # A sheet lying on a solid's face (a patch on its substrate) ties with that face: the sheet wins.
        sheets = vtk.vtkCellPicker()
        sheets.SetTolerance(0.0005)
        sheets.PickFromListOn()
        for sheet in self._sheet_actors:
            sheets.AddPickList(sheet)
        if sheets.Pick(x, y, 0, self.renderer):
            gap = np.linalg.norm(np.subtract(sheets.GetPickPosition(), picker.GetPickPosition()))
            if gap <= 1e-5 * self.camera.distance:
                return sheets.GetActor()
        return actor

    # ------------------------------------------------------------------ rendering
    def rebuild(self, project: Project, errors: list | None = None, keep_camera=True):
        cam = self.camera_position if keep_camera and self._actors else None
        for actor in self._actors.values():
            self.remove_actor(actor, render=False)
        self.remove_actor("region", render=False)
        self.remove_bounds_axes()
        self._actors.clear()
        self._tessellated.clear()
        self._colors.clear()
        self._actor_names.clear()
        self._sheet_actors.clear()
        try:
            variables = project.resolved_variables()
        except Exception as exc:  # noqa: BLE001
            if errors is not None:
                errors.append(str(exc))
            self.render()
            return

        for obj in project.objects:
            if not obj.visible:
                continue
            try:
                mesh = primitive_mesh(obj, variables, project)
            except Exception as exc:  # noqa: BLE001
                if errors is not None:
                    errors.append(f"{obj.name}: {exc}")
                continue
            mat = project.material(obj.material)
            actor = self.add_mesh(mesh, color=mat.color, opacity=1.0 - obj.transparency,
                                  show_edges=False, smooth_shading=obj.kind in ("cylinder", "sphere", "cad"),
                                  specular=0.4 if mat.conductor else 0.1, render=False,
                                  name=f"obj::{obj.name}", pickable=True)
            self._actors[obj.name] = actor
            self._tessellated[obj.name] = obj.kind == "cad"
            self._colors[obj.name] = mat.color
            # Imported (cad) objects carry no sheet flag; a flat bounding box identifies an imported sheet.
            ext = np.subtract(mesh.bounds[1::2], mesh.bounds[0::2])
            if obj.kind == "rectangle" or (obj.kind == "cad" and ext.min() <= 1e-9 * max(ext.max(), 1e-30)):
                self._sheet_actors.append(actor)
            feature_edges = mesh.extract_feature_edges(45)
            if feature_edges.n_points:  # spheres have no sharp edges
                self._actors[f"edge::{obj.name}"] = self.add_mesh(
                    feature_edges, color="#202020", line_width=1, render=False,
                    name=f"edge::{obj.name}", pickable=False)
            self._actor_names[id(actor)] = obj.name

        from .solver import model_bounds, region_bounds
        mlo, mhi = model_bounds(project, variables)
        size = float(np.linalg.norm(mhi - mlo)) or 1.0
        for port in project.ports:
            try:
                mesh = port_mesh(port, variables, size)
            except Exception as exc:  # noqa: BLE001
                if errors is not None:
                    errors.append(f"{port.name}: {exc}")
                continue
            actor = self.add_mesh(mesh, color=PORT_COLOR, opacity=1.0, render=False,
                                  name=f"port::{port.name}", pickable=True)
            self._actors[port.name] = actor
            self._actor_names[id(actor)] = port.name

        # simulation region (wireframe), like the air box in commercial tools
        try:
            lo, hi = region_bounds(project, variables)
            region = pv.Box(bounds=(lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]))
            self.add_mesh(region, style="wireframe", color="#5070a0", line_width=1, opacity=0.6,
                          name="region", pickable=False, render=False)
        except Exception:  # noqa: BLE001
            pass
        if self._show_grid and self._actors:
            self.show_bounds(bounds=(mlo[0], mhi[0], mlo[1], mhi[1], mlo[2], mhi[2]), grid="back", location="outer", color="#606060", font_size=9,
                             xtitle=f"X ({project.units})", ytitle=f"Y ({project.units})",
                             ztitle=f"Z ({project.units})", render=False)

        self.highlight(self._selected)
        if cam is not None:
            self.camera_position = cam
        else:
            self.view_isometric()
            self.fit_all()
        self.render()

    def highlight(self, name: str | None):
        self._selected = name
        for key, actor in self._actors.items():
            if key.startswith("edge::"):
                continue
            prop = actor.GetProperty()
            if self._tessellated.get(key):
                # imported geometry is triangulated: tint instead of drawing every triangle edge
                base = self._colors.get(key, (0.6, 0.6, 0.6))
                sel = pv.Color(SELECT_COLOR).float_rgb
                prop.SetColor(*(tuple(0.45 * b + 0.55 * s for b, s in zip(base, sel)) if key == name else base))
                continue
            if key == name:
                prop.SetEdgeVisibility(True)
                prop.SetEdgeColor(pv.Color(SELECT_COLOR).float_rgb)
                prop.SetLineWidth(3)
                prop.SetAmbient(0.35)
            else:
                prop.SetEdgeVisibility(False)
                prop.SetAmbient(0.0)
        self.render()

    def set_grid_visible(self, visible: bool):
        self._show_grid = visible
        if not visible:
            self.remove_bounds_axes()
        self.render()

    def fit_all(self):
        """Zoom to the model (not the much larger simulation region)."""
        actors = [a for k, a in self._actors.items() if not k.startswith("edge::")]
        if actors:
            b = np.array([a.GetBounds() for a in actors])
            self.reset_camera(bounds=(b[:, 0].min(), b[:, 1].max(), b[:, 2].min(),
                                      b[:, 3].max(), b[:, 4].min(), b[:, 5].max()))
        else:
            self.reset_camera()
        self.render()
