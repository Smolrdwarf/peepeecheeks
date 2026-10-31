"""Geometry for 'cad' primitives: shapes replayed from the project's source .aedt design.

Replay results and tessellations are cached and keyed on the source file,
its modification time and the current design-variable values, so editing a
variable in FieldForge re-runs the HFSS modeler history with the new value.
"""
from __future__ import annotations

import os

import numpy as np

_cache: dict = {"key": None, "result": None, "error": None, "meshes": {}}


def _key(project):
    src = project.source
    try:
        mtime = os.path.getmtime(src["path"])
    except OSError:
        mtime = None
    return (src.get("path"), src.get("design"), mtime, tuple(sorted(project.resolved_variables().items())))


def replay(project):
    """ReplayResult for the project's source design (cached). Raises if the source is unavailable."""
    if not project.source:
        raise RuntimeError("project has no imported CAD source")
    path = project.source.get("path", "")
    if not os.path.exists(path):
        raise FileNotFoundError(f"source design '{path}' not found (imported geometry needs the .aedt file)")
    key = _key(project)
    if _cache["key"] != key:
        from .aedt.importer import replay_for_project
        _cache.update(key=key, result=None, error=None, meshes={})
        try:
            _cache["result"] = replay_for_project(project)
        except Exception as exc:  # noqa: BLE001
            _cache["error"] = exc
    if _cache["error"] is not None:
        raise _cache["error"]
    return _cache["result"]


def shape(project, obj):
    res = replay(project)
    for pid, s in res.shapes.items():
        if res.parts[pid].name == obj.part:
            return s
    raise KeyError(f"part '{obj.part}' is not present in the source design")


def bounds(project, obj):
    if project is None:
        raise RuntimeError("bounds of a CAD object need the project")
    from .aedt.geometry import bounds as shape_bounds
    return shape_bounds(shape(project, obj))


def kind(project, obj) -> str:
    from .aedt.geometry import shape_kind
    return shape_kind(shape(project, obj))


def mesh(project, obj, slab: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Triangle mesh of an object. With slab > 0, sheets are given that thickness (closed polyhedra)."""
    from .aedt import geometry as G
    replay(project)
    key = (obj.part, round(slab, 12))
    hit = _cache["meshes"].get(key)
    if hit is not None:
        return hit
    s = shape(project, obj)
    if slab > 0 and G.shape_kind(s) == "sheet":
        s = G.sheet_to_slab(s, slab)
    lo, hi = G.bounds(s)
    deflection = max(float(np.linalg.norm(hi - lo)) * 2e-3, 1e-6)
    out = G.tessellate(s, deflection)
    _cache["meshes"][key] = out
    return out


def vertex_coords(project, obj, limit: int = 64) -> list[np.ndarray]:
    """Distinct x/y/z coordinates of an object's B-rep vertices (for mesh-line seeding)."""
    from OCP.BRep import BRep_Tool
    from OCP.TopAbs import TopAbs_VERTEX
    from OCP.TopoDS import TopoDS
    from .aedt.geometry import _iter
    pts = []
    for v in _iter(shape(project, obj), TopAbs_VERTEX):
        p = BRep_Tool.Pnt_s(TopoDS.Vertex(v))
        pts.append((p.X(), p.Y(), p.Z()))
    if not pts:
        return [np.array([])] * 3
    arr = np.unique(np.round(np.array(pts), 9), axis=0)
    out = []
    for i in range(3):
        c = np.unique(arr[:, i])
        out.append(c if len(c) <= limit else np.array([]))
    return out
