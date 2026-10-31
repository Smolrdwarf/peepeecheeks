"""Convert an Ansys HFSS design (.aedt/.aedtz) into a FieldForge Project.

What is translated, and how:

* design and project variables  -> FieldForge design variables (in model units), still editable
* parametric geometry history   -> replayed with OpenCascade (see geometry.py); each final
                                   object becomes a 'cad' primitive that is re-built from the
                                   .aedt whenever a variable changes
* materials                     -> permittivity, permeability, loss tangent; good conductors
                                   (sigma > 1e5 S/m or SolveInside=false) become PEC
* Radiation boundary on region  -> simulation region with PML on all faces
* Perfect E / Finite Cond.      -> PEC (conductor loss is not modelled)
* Lumped Port                   -> lumped port spanning the port sheet along its integration line
* Wave Port                     -> approximated by a lumped port on the port face (flagged)
* Driven setup + sweep          -> FieldForge setup (sweep range, solution frequency for far field)
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, field

import numpy as np

from ..model import (
    UNITS, ExpressionError, LumpedPort, Material, Primitive, Project, Region, Setup, evaluate,
)
from . import geometry as G
from .parser import Block, design_blocks, floats, ints, parse, unquote
from .units import convert_expr, split_value_unit, var_name

CONDUCTOR_SIGMA = 1e5     # S/m; above this a material is treated as a perfect conductor


@dataclass
class ImportReport:
    source: str
    design: str
    solution_type: str = ""
    objects: list[str] = field(default_factory=list)
    skipped_objects: list[str] = field(default_factory=list)
    variables: int = 0
    materials: list[str] = field(default_factory=list)
    ports: list[str] = field(default_factory=list)
    boundaries: list[str] = field(default_factory=list)
    setup: str = ""
    notes: list[str] = field(default_factory=list)          # translation notes (approximations)
    issues: list = field(default_factory=list)              # geometry.Issue
    op_counts: dict = field(default_factory=dict)

    @property
    def errors(self):
        return [i for i in self.issues if i.level == "error"]

    @property
    def warnings(self):
        return [i for i in self.issues if i.level != "error"]

    def summary(self) -> str:
        n_ops = sum(self.op_counts.values())
        return (f"Imported '{self.design}' from {os.path.basename(self.source)}: "
                f"{len(self.objects)} objects ({n_ops} modeler operations replayed), "
                f"{self.variables} variables, {len(self.ports)} port(s); "
                f"{len(self.errors)} error(s), {len(self.warnings) + len(self.notes)} note(s)")

    def to_html(self) -> str:
        def ul(items):
            return "<ul>" + "".join(f"<li>{_esc(i)}</li>" for i in items) + "</ul>" if items else "<p><i>none</i></p>"
        h = [f"<h2>Import report: {_esc(self.design)}</h2>",
             f"<p><b>Source:</b> {_esc(self.source)}<br><b>HFSS solution type:</b> {_esc(self.solution_type)}<br>"
             f"<b>Modeler operations replayed:</b> {sum(self.op_counts.values())}</p>"]
        if self.errors:
            h.append("<h3 style='color:#b00000'>Errors (geometry that could not be reproduced)</h3>")
            h.append(ul(str(i) for i in self.errors))
        h.append("<h3>Approximations and notes</h3>")
        h.append(ul(self.notes + [str(i) for i in self.warnings]))
        h.append(f"<h3>Objects ({len(self.objects)})</h3>" + ul(self.objects))
        if self.skipped_objects:
            h.append("<h3>Objects not simulated</h3>" + ul(self.skipped_objects))
        h.append("<h3>Excitations</h3>" + ul(self.ports))
        h.append("<h3>Boundaries</h3>" + ul(self.boundaries))
        h.append("<h3>Materials</h3>" + ul(self.materials))
        h.append(f"<h3>Analysis</h3><p>{_esc(self.setup)}</p>")
        return "\n".join(h)


def _esc(s) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


# --------------------------------------------------------------------------- public API

def list_designs(path: str) -> list[tuple[str, str]]:
    """[(design name, solution type)] for every HFSS design in the file."""
    proj = _parse_cached(path)
    return [(d.get("Name", "?"), d.get("SolutionType", "")) for d in design_blocks(proj)]


def import_aedt(path: str, design: str | None = None) -> tuple[Project, ImportReport]:
    proj_block = _parse_cached(path)
    designs = design_blocks(proj_block)
    if not designs:
        raise ValueError("The file contains no HFSS designs.")
    d = designs[0] if design is None else next((x for x in designs if x.get("Name") == design), None)
    if d is None:
        raise ValueError(f"Design '{design}' not found.")
    return _Importer(path, proj_block, d).run()


def replay_for_project(project: Project):
    """Re-run the geometry history of the project's source design with the project's variables."""
    src = project.source
    proj_block = _parse_cached(src["path"])
    d = next(x for x in design_blocks(proj_block) if x.get("Name") == src["design"])
    units = src.get("model_units", project.units)
    return G.Replayer(d, project.resolved_variables(), units, src["path"]).run()


_PARSE_CACHE: dict[str, tuple[float, Block]] = {}


def _parse_cached(path: str) -> Block:
    mtime = os.path.getmtime(path)
    hit = _PARSE_CACHE.get(path)
    if hit and hit[0] == mtime:
        return hit[1]
    blk = parse(path)
    _PARSE_CACHE[path] = (mtime, blk)
    return blk


# --------------------------------------------------------------------------- implementation

class _Importer:
    def __init__(self, path: str, proj: Block, design: Block):
        self.path = os.path.abspath(path)
        self.proj = proj
        self.d = design
        self.name = design.get("Name", "HFSSDesign")
        core = design.find("ModelSetup", "GeometryCore")
        self.units = (core.get("Units", "mm") if core is not None else "mm").replace("meter", "m")
        self.report = ImportReport(self.path, self.name, design.get("SolutionType", ""))

    # ------------------------------------------------------------------ variables
    def variables(self) -> dict[str, str]:
        raw: dict[str, str] = {}
        for blk in (self.proj.find("Properties"), self.d.find("ModelSetup", "Properties")):
            if blk is None:
                continue
            for args in blk.calls("VariableProp"):
                if len(args) < 4:
                    continue
                name, value = var_name(args[0]), args[3]
                raw[name] = _tidy(convert_expr(value, self.units))
        # drop variables that cannot be evaluated (and anything that depends on them)
        good: dict[str, str] = {}
        pending = dict(raw)
        for _ in range(len(pending) + 1):
            progressed = False
            for k, v in list(pending.items()):
                try:
                    evaluate(v, Project(variables=good).resolved_variables())
                except ExpressionError:
                    continue
                good[k] = v
                del pending[k]
                progressed = True
            if not progressed:
                break
        for k, v in pending.items():
            self.report.notes.append(f"Variable '{k}' = '{v}' could not be evaluated and was dropped")
        self.report.variables = len(good)
        return good

    # ------------------------------------------------------------------ materials
    def materials(self, used: set[str]) -> dict[str, Material]:
        out: dict[str, Material] = {}
        defs = self.proj.find("Definitions", "Materials")
        blocks = {b.name: b for b in defs.children} if defs is not None else {}
        for name in sorted(used):
            b = blocks.get(name)
            if b is None:
                self.report.notes.append(f"Material '{name}' is not defined in the project; using vacuum")
                out[name] = Material(name)
                continue
            epsr = self._mat_value(b, "permittivity", 1.0)
            mur = self._mat_value(b, "permeability", 1.0)
            tand = self._mat_value(b, "dielectric_loss_tangent", 0.0)
            sigma = self._mat_value(b, "conductivity", 0.0)
            color = (0.6, 0.6, 0.6)
            app = b.find("AttachedData", "MatAppearanceData")
            if app is not None:
                try:
                    color = tuple(int(app.get(k)) / 255 for k in ("Red", "Green", "Blue"))
                except (TypeError, ValueError):
                    pass
            conductor = sigma >= CONDUCTOR_SIGMA
            if 0 < sigma < CONDUCTOR_SIGMA:
                # bulk conductivity of a dielectric -> equivalent loss tangent is frequency dependent;
                # fold it into tan(delta) at the solution frequency.
                f = self.solution_freq_ghz * 1e9
                tand += sigma / (2 * math.pi * f * 8.854e-12 * max(epsr, 1e-9))
                self.report.notes.append(f"Material '{name}': conductivity {sigma:g} S/m folded into "
                                         f"tan(delta) at {self.solution_freq_ghz:g} GHz")
            out[name] = Material(name, epsr=epsr, mur=mur, tand=tand, conductor=conductor, color=color)
            kind = "conductor (PEC)" if conductor else f"eps_r={epsr:g}, tan_d={tand:g}"
            self.report.materials.append(f"{name}: {kind}" + (f", mu_r={mur:g}" if mur != 1 else ""))
        return out

    def _mat_value(self, b: Block, key: str, default: float) -> float:
        raw = b.get(key)
        if raw is None:
            sub = b.child(key)
            if sub is not None:
                self.report.notes.append(f"Material '{b.name}': {key} is anisotropic/frequency dependent; "
                                         f"first component used")
                comp = sub.get("component1") or sub.get("Component1")
                raw = comp
        if raw is None:
            return default
        try:
            return float(evaluate(convert_expr(raw, self.units), self.var_values))
        except Exception:  # noqa: BLE001
            self.report.notes.append(f"Material '{b.name}': could not evaluate {key}='{raw}', using {default}")
            return default

    # ------------------------------------------------------------------ setup
    def setup(self) -> Setup:
        setups = self.d.find("AnalysisSetup", "SolveSetups")
        s = Setup()
        drv = None
        if setups is not None:
            drv = next((c for c in setups.children
                        if c.get("SetupType", "").startswith("HfssDriven") and c.get("IsEnabled", "true") != "false"),
                       None) or next((c for c in setups.children if c.get("SetupType", "").startswith("HfssDriven")),
                                     None)
        if drv is None:
            self.report.notes.append("No driven solution setup found; using default 1-3 GHz sweep")
            self.solution_freq_ghz = 2.0
            return s
        s.name = drv.name
        f0 = self._freq(drv.get("Frequency")) or 2.0
        self.solution_freq_ghz = f0
        lo, hi, npts = None, None, 0
        sweeps = drv.child("Sweeps")
        for sw in (sweeps.children if sweeps is not None else []):
            if sw.get("IsEnabled", "true") == "false":
                continue
            rt = sw.get("RangeType", "")
            a, b = self._freq(sw.get("RangeStart")), self._freq(sw.get("RangeEnd"))
            if rt == "SinglePoints":
                a = b = a or self._freq(sw.get("RangeStart"))
            if a is None or b is None:
                continue
            lo = a if lo is None else min(lo, a)
            hi = b if hi is None else max(hi, b)
            if rt == "LinearCount":
                npts = max(npts, int(float(sw.get("RangeCount", "401"))))
            elif rt == "LinearStep":
                step = self._freq(sw.get("RangeStep")) or (b - a) / 400
                npts = max(npts, int(round((b - a) / step)) + 1)
        if lo is None:
            lo, hi = 0.5 * f0, 1.5 * f0
            self.report.notes.append("Setup has no frequency sweep; using 0.5-1.5x the solution frequency")
        if lo <= 0:
            lo = hi / 100
        s.f_start, s.f_stop = _fmt(lo), _fmt(hi)
        s.n_points = max(npts, 101) if npts else 401
        s.farfield_freq = _fmt(min(max(f0, lo), hi))
        s.farfield = True
        # HFSS adaptive meshing -> FieldForge adaptive passes (global FDTD re-mesh, same delta-S criterion).
        # PercentRefinement (tetrahedra added per pass) has no FDTD equivalent; refine_factor keeps its default.
        def num(key, conv, default):
            try:
                return conv(float(drv.get(key)))
            except (TypeError, ValueError):
                return default
        s.adaptive = True
        s.max_passes = num("MaximumPasses", int, s.max_passes)
        s.max_delta_s = num("MaxDeltaS", float, s.max_delta_s)
        s.min_converged_passes = num("MinimumConvergedPasses", int, s.min_converged_passes)
        self.report.setup = (f"{drv.name}: sweep {lo:g}-{hi:g} GHz ({s.n_points} points), solution frequency "
                             f"{f0:g} GHz (used for the far field). Adaptive meshing: up to {s.max_passes} passes "
                             f"until max delta S < {s.max_delta_s:g} ({s.min_converged_passes} converged pass(es)), "
                             f"as in the HFSS setup; FieldForge re-meshes the whole FDTD grid "
                             f"{s.refine_factor:g}x finer per pass, starting at {s.cells_per_wavelength} cells per "
                             f"wavelength.")
        return s

    def _freq(self, raw) -> float | None:
        if not raw:
            return None
        try:
            return float(evaluate(convert_expr(raw, self.units), self.var_values))
        except Exception:  # noqa: BLE001
            return None

    # ------------------------------------------------------------------ main
    def run(self) -> tuple[Project, ImportReport]:
        r = self.report
        variables = self.variables()
        self.var_values = Project(variables=variables).resolved_variables()
        setup = self.setup()

        res = G.Replayer(self.d, self.var_values, self.units, self.path).run()
        r.issues = list(res.issues)
        r.op_counts = res.op_counts
        if res.unsupported:
            r.notes.append("Unsupported modeler operations: " +
                           ", ".join(f"{k} (x{v})" for k, v in sorted(res.unsupported.items())))

        units = self.units if self.units in UNITS else "mm"
        if units != self.units:
            r.notes.append(f"Model units '{self.units}' shown as mm")
        p = Project(name=_safe_name(self.name), units=units, variables=variables, setup=setup)
        p.source = {"type": "aedt", "path": self.path, "design": self.name, "model_units": self.units}

        bounds = self.d.find("BoundarySetup", "Boundaries")
        bnds = list(bounds.children) if bounds is not None else []
        port_parts: set[int] = set()
        pec_parts: set[int] = set()
        radiation_parts: set[int] = set()
        for b in bnds:
            bt = b.get("BoundType", "")
            objs = ints(b.call("Objects") or [])
            if bt in ("Lumped Port", "Wave Port"):
                port_parts.update(objs)
            elif bt in ("Perfect E", "Finite Conductivity", "Layered Impedance", "Impedance"):
                pec_parts.update(objs)
            elif bt == "Radiation":
                radiation_parts.update(objs)

        # materials actually used by simulated objects
        used = set()
        for pid, info in res.parts.items():
            if pid in res.shapes and info.material:
                used.add(info.material)
        self.solution_freq_ghz = getattr(self, "solution_freq_ghz", 2.0)
        mats = self.materials(used)
        p.materials = mats

        # objects
        for pid, shape in res.shapes.items():
            info = res.parts[pid]
            kind = G.shape_kind(shape)
            label = f"{info.name} ({kind}, {info.material or 'no material'})"
            if pid in port_parts:
                r.skipped_objects.append(f"{info.name}: port sheet (represented by the port)")
                continue
            if pid in radiation_parts:
                r.skipped_objects.append(f"{info.name}: radiation boundary box (used as the simulation region)")
                continue
            if not info.model:
                r.skipped_objects.append(f"{info.name}: non-model object")
                continue
            if kind == "wire":
                r.skipped_objects.append(f"{info.name}: line object (no EM effect)")
                continue
            mat = info.material or "vacuum"
            if pid in pec_parts:
                mat = "pec"
            elif mat in mats and mats[mat].conductor:
                pass
            elif not info.solve_inside and mat in mats:
                # HFSS treats SolveInside=false objects as conductors (field excluded)
                mats[mat].conductor = True
            elif kind == "sheet":
                # a sheet with no boundary assigned is electromagnetically transparent in HFSS
                r.skipped_objects.append(f"{info.name}: sheet without boundary (transparent in HFSS)")
                obj = Primitive(_safe_name(info.name), "cad", "vacuum", transparency=0.85)
                obj.part = info.name
                obj.visible = True
                p.objects.append(obj)
                continue
            obj = Primitive(_safe_name(info.name), "cad", mat,
                            transparency=min(max(info.transparency, 0.0), 0.9))
            obj.part = info.name
            p.objects.append(obj)
            r.objects.append(label + (" -> PEC" if mat == "pec" else ""))

        # region / boundaries
        if res.region is not None:
            lo, hi = res.region["lo"], res.region["hi"]
            p.region = Region(padding="0", boundaries=["PML_8"] * 6)
            p.region.explicit = [float(lo[0]), float(hi[0]), float(lo[1]), float(hi[1]), float(lo[2]), float(hi[2])]
        for b in bnds:
            bt = b.get("BoundType", "")
            if bt == "Radiation":
                r.boundaries.append(f"{b.name}: Radiation -> PML (8 cells) on the region faces")
                if res.region is None:
                    for pid in ints(b.call("Objects") or []):
                        if pid in res.shapes:
                            lo, hi = G.bounds(res.shapes[pid])
                            p.region.explicit = [float(lo[0]), float(hi[0]), float(lo[1]), float(hi[1]),
                                                 float(lo[2]), float(hi[2])]
            elif bt in ("Perfect E", "Finite Conductivity"):
                extra = " (conductor loss neglected)" if bt == "Finite Conductivity" else ""
                if b.call("Faces"):
                    r.notes.append(f"Boundary '{b.name}' ({bt}) is assigned to faces; only object assignments "
                                   "are translated")
                r.boundaries.append(f"{b.name}: {bt} -> PEC{extra}")
            elif bt in ("Lumped Port", "Wave Port", "Terminal"):
                continue
            else:
                r.notes.append(f"Boundary '{b.name}' of type '{bt}' is not supported and was ignored")

        # ports
        self.conductor_parts = set(pec_parts) | {
            pid for pid, info in res.parts.items()
            if not info.solve_inside or (info.material in mats and mats[info.material].conductor)}
        self._ports(p, res, bnds)
        if not p.ports:
            r.notes.append("No supported excitation was found; add a lumped port before simulating")
        try:
            from ..solver import port_gaps
            r.notes += port_gaps(p)
        except Exception as exc:  # noqa: BLE001  (never fail an import over a diagnostic)
            r.notes.append(f"Port connectivity check skipped: {exc}")
        return p, r

    # ------------------------------------------------------------------ ports
    def _ports(self, p: Project, res, bnds):
        r = self.report
        first = True
        for b in bnds:
            bt = b.get("BoundType", "")
            if bt not in ("Lumped Port", "Wave Port"):
                continue
            z = 50.0
            if b.get("Impedance"):
                try:
                    z = float(evaluate(convert_expr(b.get("Impedance"), self.units), self.var_values))
                except Exception:  # noqa: BLE001
                    pass
            if b.get("FullResistance"):
                try:
                    z = float(evaluate(convert_expr(b.get("FullResistance"), self.units), self.var_values))
                except Exception:  # noqa: BLE001
                    pass
            lo, hi = self._port_extent(b, res)
            if lo is None:
                r.notes.append(f"Port '{b.name}' ({bt}): could not determine its location; skipped")
                continue
            int_line = self._int_line(b)
            axis, (lo, hi), approx = self._port_box(lo, hi, int_line, res, bt)
            port = LumpedPort(_safe_name(b.name), axis=axis, impedance=_fmt(z), excited=first,
                              params={"x1": _fmt(lo[0]), "y1": _fmt(lo[1]), "z1": _fmt(lo[2]),
                                      "x2": _fmt(hi[0]), "y2": _fmt(hi[1]), "z2": _fmt(hi[2])})
            p.ports.append(port)
            first = False
            desc = f"{b.name}: {bt} -> lumped port, {z:g} Ohm, current along {axis}"
            if approx:
                desc += f" ({approx})"
                r.notes.append(f"Port '{b.name}': {approx}")
            r.ports.append(desc)
        if len(p.ports) > 1:
            r.notes.append("FieldForge excites one port per run; only the first port is excited, "
                           "giving the first column of the S-matrix")

    def _port_extent(self, b: Block, res):
        objs = ints(b.call("Objects") or [])
        for pid in objs:
            s = res.shapes.get(pid)
            if s is not None:
                return G.bounds(s)
        faces = ints(b.call("Faces") or [])
        rep = G.Replayer(self.d, self.var_values, self.units)
        for fid in faces:
            ref = rep.face_ref(fid)
            if ref and ref.get("verts"):
                v = np.array(ref["verts"])
                return v.min(axis=0), v.max(axis=0)
            if ref and ref.get("mid"):
                # locate the face on the final shapes
                mid = G.vertex(ref["mid"])
                for s in res.shapes.values():
                    for f in G.faces_of(s):
                        if G.dist(mid, f) < 1e-6:
                            return G.bounds(f)
        return None, None

    def _int_line(self, b: Block):
        mode = b.find("Modes", "Mode1")
        il = mode.child("IntLine") if mode is not None else None
        if il is None or mode.get("UseIntLine") == "false":
            return None
        pts = []
        for gp in il.children_named("GeometryPosition"):
            try:
                pts.append([self._len(gp.get(k)) for k in ("XPosition", "YPosition", "ZPosition")])
            except Exception:  # noqa: BLE001
                return None
        return pts if len(pts) == 2 else None

    def _len(self, raw) -> float:
        return float(evaluate(convert_expr(raw, self.units), self.var_values))

    def _port_box(self, lo, hi, int_line, res, bt):
        lo, hi = np.array(lo, float), np.array(hi, float)
        ext = hi - lo
        size = max(ext.max(), 1e-9)
        approx = ""
        if int_line is not None:
            a, b = np.array(int_line[0]), np.array(int_line[1])
            # AEDT occasionally stores one end in a different scale; snap it back onto the sheet
            for pt in (a, b):
                for _ in range(2):
                    if np.any(pt < lo - size) or np.any(pt > hi + size):
                        pt /= 1000.0
            d = b - a
            axis = int(np.argmax(np.abs(d)))
            if np.linalg.norm(d) > 0 and abs(d[axis]) / np.linalg.norm(d) < 0.98:
                approx = "integration line is not axis aligned; snapped to the nearest axis"
            lo2, hi2 = lo.copy(), hi.copy()
            lo2[axis], hi2[axis] = sorted((a[axis], b[axis]))
            if bt == "Wave Port":
                flat = int(np.argmin(ext))
                narrowed, _ = self._shrink_wave_port(lo2, hi2, flat, axis, res, keep_axis=True), None
                lo2, hi2 = narrowed
                approx = ("wave port approximated by a lumped port along its integration line, narrowed to the "
                          "signal conductor; expect differences from HFSS modal port de-embedding")
            return "XYZ"[axis], (lo2, hi2), approx
        # no integration line: current flows along the face's in-plane axis that crosses
        # from the ground to the signal conductor; guess the vertical (Z) one if present
        flat = int(np.argmin(ext))
        in_plane = [i for i in range(3) if i != flat]
        axis = 2 if 2 in in_plane else in_plane[int(np.argmin([ext[i] for i in in_plane]))]
        lo2, hi2 = lo.copy(), hi.copy()
        if bt == "Wave Port":
            lo2, hi2 = self._shrink_wave_port(lo, hi, flat, axis, res)
            approx = ("wave port approximated by a lumped port between the conductors touching the port face; "
                      "expect differences from HFSS modal port de-embedding")
        else:
            approx = "no integration line stored; current direction guessed"
        return "XYZ"[axis], (lo2, hi2), approx

    def _shrink_wave_port(self, lo, hi, flat, axis, res, keep_axis=False):
        """Shrink a wave-port face to the gap between ground and signal conductor touching it."""
        plane = 0.5 * (lo[flat] + hi[flat])
        tol = 1e-6 + 1e-6 * abs(plane)
        touching = []
        for pid, s in res.shapes.items():
            if pid not in getattr(self, "conductor_parts", set()):
                continue
            a, b = G.bounds(s)
            if a[flat] - tol <= plane <= b[flat] + tol and np.all(b >= lo - tol) and np.all(a <= hi + tol):
                touching.append((a, b))
        lo2, hi2 = lo.copy(), hi.copy()
        if len(touching) >= 2 or (keep_axis and touching):
            touching.sort(key=lambda ab: ab[0][axis])
            ground, signal = touching[0], touching[-1]
            if not keep_axis:
                lo2[axis] = ground[1][axis]
                hi2[axis] = signal[0][axis]
            other = [i for i in range(3) if i not in (flat, axis)][0]
            lo2[other] = max(lo[other], signal[0][other])
            hi2[other] = min(hi[other], signal[1][other])
            if hi2[axis] <= lo2[axis] or hi2[other] <= lo2[other]:
                return lo, hi
        return lo2, hi2


def _tidy(expr: str) -> str:
    """'(60.0)' -> '60' for readability."""
    s = expr.strip()
    while s.startswith("(") and s.endswith(")") and s.count("(") == 1:
        s = s[1:-1].strip()
    try:
        v = float(s)
        return _fmt(v)
    except ValueError:
        return s


def _fmt(v: float) -> str:
    return f"{float(v):.10g}"


def _safe_name(name: str) -> str:
    out = "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in name)
    return out if out and not out[0].isdigit() else "_" + out
