"""Replay the parametric modeler history of an AEDT design with OpenCascade.

AEDT stores geometry as a time-ordered list of operations (create box, cover
lines, move, subtract, ...) spread across the parts they act on. Operations
from different parts interleave (a Clone copies its source *as it was at that
moment*), so every operation is gathered, sorted by its global ID and executed
in order against a ``part id -> TopoDS_Shape`` table.

Faces and edges are referenced by kernel IDs. AEDT also stores the geometry
of referenced entities (vertices, mid point, area), which is used to find the
matching face/edge on the replayed shape, and to rebuild planar faces of
imported (Parasolid) bodies that cannot be read directly.
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, field

import numpy as np
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_Copy, BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakePolygon,
    BRepBuilderAPI_MakeWire, BRepBuilderAPI_RightCorner, BRepBuilderAPI_Transform,
)
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepFilletAPI import BRepFilletAPI_MakeChamfer, BRepFilletAPI_MakeFillet
from OCP.BRepGProp import BRepGProp
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRepOffsetAPI import BRepOffsetAPI_MakePipeShell
from OCP.BRepPrimAPI import (
    BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone, BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakePrism,
    BRepPrimAPI_MakeSphere,
)
from OCP.Bnd import Bnd_Box
from OCP.BRepLib import BRepLib_FindSurface
from OCP.GC import GC_MakeArcOfCircle, GC_MakeCircle
from OCP.GeomAPI import GeomAPI_Interpolate
from OCP.GProp import GProp_GProps
from OCP.STEPControl import STEPControl_Reader
from OCP.collections import HArray1_gp_Pnt as TColgp_HArray1OfPnt
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_REVERSED, TopAbs_SHELL, TopAbs_SOLID, TopAbs_VERTEX, TopAbs_WIRE
from OCP.TopExp import TopExp_Explorer
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS, TopoDS_Compound, TopoDS_Shape
from OCP.gp import gp_Ax1, gp_Ax2, gp_Ax3, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec

from ..model import ExpressionError, evaluate
from .parser import Block, floats, ints, split_args
from .units import convert_expr

AXIS_VEC = {"X": (1, 0, 0), "Y": (0, 1, 0), "Z": (0, 0, 1)}
# For a sheet with the given normal: (first in-plane axis, second in-plane axis)
PLANE_AXES = {"Z": ("X", "Y"), "X": ("Y", "Z"), "Y": ("Z", "X")}


# --------------------------------------------------------------------------- data

@dataclass
class PartInfo:
    part_id: int
    name: str
    material: str = ""
    solve_inside: bool = True
    color: tuple = (0.6, 0.6, 0.6)
    transparency: float = 0.0
    model: bool = True          # False for 'NonModel' objects (construction geometry)
    toplevel: bool = False
    group: str = ""


@dataclass
class Issue:
    level: str                  # 'warn' | 'error'
    part: str
    op: str
    message: str

    def __str__(self):
        where = f"{self.part}: " if self.part else ""
        return f"{where}{self.op}: {self.message}" if self.op else f"{where}{self.message}"


@dataclass
class Op:
    type: str
    id: int
    part_id: int
    cs: int
    block: Block
    params: Block | None
    owner: int = -1             # part whose operation list holds this op
    seq: int = 0                # position within that list


@dataclass
class ReplayResult:
    shapes: dict[int, TopoDS_Shape]
    parts: dict[int, PartInfo]
    region: dict | None
    issues: list[Issue] = field(default_factory=list)
    op_counts: dict[str, int] = field(default_factory=dict)
    unsupported: dict[str, int] = field(default_factory=dict)


# --------------------------------------------------------------------------- helpers

def _iter(shape, kind):
    exp = TopExp_Explorer(shape, kind)
    while exp.More():
        yield exp.Current()
        exp.Next()


def faces_of(shape):
    return [TopoDS.Face(f) for f in _iter(shape, TopAbs_FACE)]


def edges_of(shape):
    return [TopoDS.Edge(e) for e in _iter(shape, TopAbs_EDGE)]


def count(shape, kind) -> int:
    return sum(1 for _ in _iter(shape, kind))


def is_null(shape) -> bool:
    return shape is None or shape.IsNull()


def compound(shapes) -> TopoDS_Shape | None:
    shapes = [s for s in shapes if not is_null(s)]
    if not shapes:
        return None
    if len(shapes) == 1:
        return shapes[0]
    comp = TopoDS_Compound()
    b = BRep_Builder()
    b.MakeCompound(comp)
    for s in shapes:
        b.Add(comp, s)
    return comp


def transformed(shape, trsf: gp_Trsf, copy=True):
    return BRepBuilderAPI_Transform(shape, trsf, copy).Shape()


def bounds(shape) -> tuple[np.ndarray, np.ndarray]:
    box = Bnd_Box()
    # exact box: Add_s pads it by the shape tolerance (1e-7), which puts mesh lines a hair inside
    # neighbouring conductors and can disable a lumped port
    BRepBndLib.AddOptimal_s(shape, box, False, False)
    if box.IsVoid():
        return np.zeros(3), np.zeros(3)
    lo, hi = box.CornerMin(), box.CornerMax()
    return np.array([lo.X(), lo.Y(), lo.Z()]), np.array([hi.X(), hi.Y(), hi.Z()])


def face_area(face) -> float:
    props = GProp_GProps()
    BRepGProp.SurfaceProperties_s(face, props)
    return props.Mass()


def face_normal(face) -> np.ndarray | None:
    """Normal of a planar face (respecting face orientation)."""
    surf = BRepLib_FindSurface(face, 1e-6, True)
    if not surf.Found():
        return None
    from OCP.GeomLProp import GeomLProp_SLProps
    s = BRep_Tool.Surface_s(face)
    u0, u1, v0, v1 = _uv_bounds(face)
    props = GeomLProp_SLProps(s, 0.5 * (u0 + u1), 0.5 * (v0 + v1), 1, 1e-6)
    if not props.IsNormalDefined():
        return None
    n = props.Normal()
    v = np.array([n.X(), n.Y(), n.Z()])
    if face.Orientation() == TopAbs_REVERSED:
        v = -v
    return v


def _uv_bounds(face):
    from OCP.BRepTools import BRepTools
    return BRepTools.UVBounds_s(face)


def dist(shape_a, shape_b) -> float:
    d = BRepExtrema_DistShapeShape(shape_a, shape_b)
    return d.Value() if d.IsDone() else float("inf")


def vertex(p) -> TopoDS_Shape:
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    return BRepBuilderAPI_MakeVertex(gp_Pnt(*p)).Vertex()


def _pnt(p):
    return gp_Pnt(float(p[0]), float(p[1]), float(p[2]))


# --------------------------------------------------------------------------- replay engine

class Replayer:
    def __init__(self, design: Block, variables: dict[str, float], units: str, aedt_path: str = "",
                 step_search_dirs: list[str] | None = None):
        self.design = design
        self.vars = variables
        self.units = units
        self.aedt_path = aedt_path
        self.step_dirs = step_search_dirs or []
        self.shapes: dict[int, TopoDS_Shape] = {}
        self.parts: dict[int, PartInfo] = {}
        self.cs: dict[int, gp_Trsf] = {1: gp_Trsf()}   # local -> global
        self.issues: list[Issue] = []
        self.pending: dict[str, list] = {}               # To-op family -> produced shapes
        self.region_params: Block | None = None
        self.region_part: int | None = None
        self.op_counts: dict[str, int] = {}
        self.unsupported: dict[str, int] = {}
        self._face_index: dict[int, dict] | None = None
        self._edge_index: dict[int, dict] | None = None
        self._step_cache: dict[str, list] = {}
        # Rigid transforms applied to each part since creation. Stored face/edge
        # reference geometry is recorded at creation time, so it must be moved too.
        self.part_trsf: dict[int, gp_Trsf] = {}

    # ------------------------------------------------------------------ expressions
    def num(self, expr, default=None) -> float:
        if expr is None:
            if default is None:
                raise ExpressionError("missing value")
            return default
        return evaluate(convert_expr(expr, self.units), self.vars)

    def p(self, params: Block, key: str, default=None) -> float:
        return self.num(params.get(key), default)

    def warn(self, part_id, op, msg, level="warn"):
        name = self.parts[part_id].name if part_id in self.parts else ""
        self.issues.append(Issue(level, name, op, msg))

    # ------------------------------------------------------------------ coordinate systems
    def to_global(self, shape, cs_id):
        t = self.cs.get(cs_id)
        if t is None or cs_id == 1 or shape is None:
            return shape
        return transformed(shape, t)

    def cs_point(self, p, cs_id):
        t = self.cs.get(cs_id, gp_Trsf())
        q = _pnt(p).Transformed(t)
        return np.array([q.X(), q.Y(), q.Z()])

    def cs_vec(self, v, cs_id):
        t = self.cs.get(cs_id, gp_Trsf())
        w = gp_Vec(*map(float, v)).Transformed(t)
        return np.array([w.X(), w.Y(), w.Z()])

    def _make_cs(self, op: Op):
        prm = op.params
        mode = prm.get("Mode", "Axis/Position")
        if mode != "Axis/Position":
            self.warn(-1, "CoordinateSystem", f"mode '{mode}' not supported; using reference CS")
            self.cs[op.id] = self.cs.get(op.cs, gp_Trsf())
            return
        o = [self.p(prm, k) for k in ("OriginX", "OriginY", "OriginZ")]
        x = [self.p(prm, k) for k in ("XAxisXvec", "XAxisYvec", "XAxisZvec")]
        y = [self.p(prm, k) for k in ("YAxisXvec", "YAxisYvec", "YAxisZvec")]
        z = np.cross(x, y)
        ax3 = gp_Ax3(_pnt(o), gp_Dir(*z), gp_Dir(*x))
        local = gp_Trsf()
        local.SetTransformation(ax3)      # global -> local
        local.Invert()                    # local -> global (relative to reference CS)
        parent = self.cs.get(op.cs, gp_Trsf())
        t = gp_Trsf()
        t.Multiply(parent)
        t.Multiply(local)
        self.cs[op.id] = t

    # ------------------------------------------------------------------ entity indices
    def _build_indices(self):
        faces, edges = {}, {}

        def collect(b: Block, parent_kid):
            kid = b.props.get("kID") or b.props.get("ID") or parent_kid
            fg = b.child("FaceGeomTopol")
            if fg is not None and kid is not None:
                geo = fg.child("FaceGeometry")
                if geo is not None:
                    mid = geo.call("FcUVMid")
                    verts = [floats(v)[:3] for v in (geo.child("FcTolVts").calls("TolVt") if geo.child("FcTolVts") else [])]
                    topo = fg.call("FaceTopol")
                    faces.setdefault(int(kid), {
                        "mid": floats(mid) if mid else None, "area": float(geo.get("Area", "0") or 0),
                        "verts": verts, "loops": int(topo[0]) if topo else 1})
            eg = b.child("EdgeGeomTopol")
            if eg is not None and kid is not None:
                mid = eg.call("EdgeMidPoint")
                verts = [floats(v)[:3] for v in (eg.child("EdTolVts").calls("TolVt") if eg.child("EdTolVts") else [])]
                edges.setdefault(int(kid), {"mid": floats(mid) if mid else None, "verts": verts})
            if b.name in ("Edge",) and kid is not None and b.call("EdgeMidPoint"):
                verts = [floats(v)[:3] for v in (b.child("EdTolVts").calls("TolVt") if b.child("EdTolVts") else [])]
                edges.setdefault(int(kid), {"mid": floats(b.call("EdgeMidPoint")), "verts": verts})
            for c in b.children:
                collect(c, b.props.get("kID", parent_kid) if b.name.startswith("Assigned") else None)

        collect(self.design, None)
        self._face_index, self._edge_index = faces, edges

    def face_ref(self, fid):
        if self._face_index is None:
            self._build_indices()
        return self._face_index.get(fid)

    def edge_ref(self, eid):
        if self._edge_index is None:
            self._build_indices()
        return self._edge_index.get(eid)

    def resolve_face(self, part_id, fid):
        """Find the face with AEDT id `fid` on the current shape of `part_id`."""
        ref = self.face_ref(fid)
        shape = self.shapes.get(part_id)
        if ref is None:
            return None
        if not is_null(shape) and ref.get("mid"):
            faces = faces_of(shape)
            areas = [face_area(f) for f in faces]
            size = np.linalg.norm(np.subtract(*bounds(shape)[::-1]))
            tol = 1e-4 * max(size, 1e-9) + 1e-7
            for mid_pt in (self._moved(part_id, ref["mid"]), ref["mid"]):
                mid = vertex(mid_pt)
                scored = [(dist(mid, f) + 1e-3 * abs(a - ref["area"]) / max(ref["area"], 1e-12), f)
                          for f, a in zip(faces, areas)]
                score, best = min(scored, key=lambda x: x[0], default=(float("inf"), None))
                if best is not None and score < tol:
                    return best
            # positions did not match (history we cannot track); accept a unique area match
            same = [f for f, a in zip(faces, areas) if abs(a - ref["area"]) <= 1e-4 * max(ref["area"], 1e-12)]
            if len(same) == 1:
                return same[0]
        # Fallback: rebuild a planar face from the stored outline (imported bodies)
        return self._face_from_outline(ref)

    def _moved(self, part_id, p):
        t = self.part_trsf.get(part_id)
        if t is None:
            return p
        q = _pnt(p).Transformed(t)
        return [q.X(), q.Y(), q.Z()]

    def _face_from_outline(self, ref):
        verts = np.array(ref.get("verts") or [])
        if len(verts) < 3 or ref.get("loops", 1) != 1:
            return None
        c = verts.mean(axis=0)
        _, s, vt = np.linalg.svd(verts - c)
        if s[-1] > 1e-6 * max(s[0], 1e-12) + 1e-9:
            return None  # not planar
        from OCP.BRepCheck import BRepCheck_Analyzer
        # AEDT lists the vertices in boundary order; only re-sort if that gives an invalid face
        for ordered in (verts, _order_loop(verts, vt[0], vt[1])):
            for with_arcs in (True, False):
                wire = _outline_wire(ordered, c, vt[0], vt[1], with_arcs)
                if wire is None:
                    continue
                f = BRepBuilderAPI_MakeFace(wire, True)
                if f.IsDone() and BRepCheck_Analyzer(f.Face()).IsValid():
                    return f.Face()
        return None

    def resolve_edge(self, part_id, eid):
        ref = self.edge_ref(eid)
        if ref is None:
            return None
        shape = self.shapes.get(part_id)
        if not is_null(shape) and ref.get("mid"):
            for mid_pt in (self._moved(part_id, ref["mid"]), ref["mid"]):
                mid = vertex(mid_pt)
                best = min(edges_of(shape), key=lambda e: dist(mid, e), default=None)
                if best is not None and dist(mid, best) < 1e-6 + 1e-6 * np.abs(mid_pt).max():
                    return best
        verts = ref.get("verts") or []
        if len(verts) == 2:
            mid = ref.get("mid")
            a, b = np.array(verts[0]), np.array(verts[1])
            if mid is not None and np.linalg.norm(np.cross(b - a, np.array(mid) - a)) > 1e-9 * (1 + np.linalg.norm(b - a)) ** 2:
                arc = GC_MakeArcOfCircle(_pnt(a), _pnt(mid), _pnt(b)).Value()
                return BRepBuilderAPI_MakeEdge(arc).Edge()
            return BRepBuilderAPI_MakeEdge(_pnt(a), _pnt(b)).Edge()
        return None

    # ------------------------------------------------------------------ gathering
    def gather(self) -> list[Op]:
        geo = self.design.find("ModelSetup", "GeometryCore", "GeometryOperations")
        if geo is None:
            raise ValueError("design has no geometry")
        ops: list[Op] = []
        csb = geo.child("CoordinateSystems")
        if csb is not None:
            for seq, b in enumerate(csb.children):
                op = self._op(b, -1)
                op.seq = seq
                ops.append(op)
        for container, top in (("ToplevelParts", True), ("OperandParts", False)):
            cb = geo.child(container)
            if cb is None:
                continue
            for part in cb.children_named("GeometryPart"):
                ob = part.child("Operations")
                if ob is None or not ob.children:
                    continue
                first = ob.children[0]
                pid = int(first.get("ParentPartID", "-1"))
                attrs = part.child("Attributes")
                self.parts[pid] = _part_info(pid, attrs, top)
                for seq, b in enumerate(ob.children):
                    op = self._op(b, pid)
                    op.seq = seq
                    ops.append(op)
        ops.sort(key=lambda o: o.id)
        return ops

    @staticmethod
    def _op(b: Block, owner) -> Op:
        params = next((c for c in b.children if c.name.endswith("Parameters")), None)
        return Op(b.get("OperationType", b.name), int(b.get("ID", "0")), int(b.get("ParentPartID", owner) or owner),
                  int(b.get("ReferenceCoordSystemID", "1") or 1), b, params, owner)

    # ------------------------------------------------------------------ main loop
    #
    # Operation IDs are only roughly chronological: duplicated geometry gets its
    # history renumbered, so e.g. a Unite can carry a smaller ID than the
    # creation of its own tool parts. Replay is therefore dependency driven:
    # each part's own list is executed in order, a boolean first completes its
    # tool parts, and a "From" operation first runs its source part up to the
    # matching "To" operation. Global ID order is only the tie-breaker.

    FROM_TO = {
        "BodyFromFaceFrom": "BodyFromFace", "BodyFromEdgeFrom": "BodyFromEdge",
        "SweepFaceAlongNormalFrom": "SweepFaceAlongNormal", "CloneFrom": "CloneTo",
        "DuplicateBodyMirror": "DuplicateMirror", "DuplicateBodyAlongLine": "DuplicateAlongLine",
        "DuplicateBodyAroundAxis": "DuplicateAroundAxis",
    }
    BOOLEANS = ("Unite", "Substract", "Subtract", "Intersect")

    def run(self) -> ReplayResult:
        ops = self.gather()
        self.by_owner: dict[int, list[Op]] = {}
        for op in ops:
            self.by_owner.setdefault(op.owner, []).append(op)
        for lst in self.by_owner.values():
            lst.sort(key=lambda o: o.seq)
        self.by_id: dict[int, list[Op]] = {}
        for op in ops:
            self.by_id.setdefault(op.id, []).append(op)
        self._done: set[int] = set()
        self._active: set[int] = set()
        self._claimed: dict[int, Op] = {}         # id(From-op) -> To-op
        for op in ops:
            self.op_counts[op.type] = self.op_counts.get(op.type, 0) + 1
        for op in ops:
            self._execute(op)
        region = None
        if self.region_params is not None:
            region = self._region_box()
        final = {pid: s for pid, s in self.shapes.items()
                 if pid in self.parts and self.parts[pid].toplevel and not is_null(s)}
        for pid, info in self.parts.items():
            if info.toplevel and pid not in final and pid != self.region_part:
                self.warn(pid, "", "no geometry could be reconstructed for this object", "error")
        # Failures on intermediate (consumed) objects only matter through their effect on final
        # objects, which is reported there; keep them as warnings.
        final_names = {self.parts[pid].name for pid in final} | {
            info.name for info in self.parts.values() if info.toplevel}
        for i in self.issues:
            if i.level == "error" and i.part and i.part not in final_names:
                i.level = "warn"
                i.message += " (intermediate object)"
        return ReplayResult(final, self.parts, region, self.issues, self.op_counts, self.unsupported)

    def _advance(self, owner: int, until: Op | None = None, inclusive=False):
        """Execute the ops of `owner` in list order, up to (optionally including) `until`."""
        for o in self.by_owner.get(owner, []):
            if o is until:
                if inclusive:
                    self._execute(o)
                return
            self._execute(o)

    def _execute(self, op: Op):
        key = id(op)
        if key in self._done or key in self._active:
            return
        self._active.add(key)
        try:
            self._advance(op.owner, until=op)
            self._resolve_dependencies(op)
            self._dispatch(op)
        finally:
            self._active.discard(key)
            self._done.add(key)

    def _resolve_dependencies(self, op: Op):
        q = op.params
        if q is None:
            return
        if op.type in self.BOOLEANS:
            for t in _idlist(q, "ToolParts"):
                self._advance(t)
        elif op.type == "SweepAlongVector":
            pid = int(q.get("ProfileID", op.part_id))
            if pid != op.owner:
                self._advance(pid)
        if op.type in self.FROM_TO:
            to_op = self._to_op_for(op)
            if to_op is not None:
                self._execute(to_op)

    def _to_op_for(self, op: Op) -> Op | None:
        if id(op) in self._claimed:
            return self._claimed[id(op)]
        family = self.FROM_TO[op.type]
        found = None
        # 1) explicit link stored in the operation identity
        for b in op.block.walk():
            for k, v in b.props.items():
                if k.endswith("ToOpnId"):
                    found = next((o for o in self.by_id.get(int(v), []) if o.type == family), None)
        q = op.params
        which = 0
        if q is not None:
            which = int(q.get("WhichClone") or q.get("WhichFace") or q.get("WhichEdge") or 0)
        # 2) source part named explicitly (clones / duplicates)
        if found is None and q is not None and q.get("SourceID"):
            src = int(q.get("SourceID"))
            used = {id(t) for t in self._claimed.values()}
            cands = [o for o in self.by_owner.get(src, []) if o.type == family]
            if which > 0:
                prev = [t for t in self._claimed.values() if t in cands]
                found = prev[-1] if prev else (cands[0] if cands else None)
            else:
                found = next((o for o in cands if id(o) not in used), cands[-1] if cands else None)
        # 3) nearest preceding To-op of the same family by ID
        if found is None:
            cands = [o for lst in self.by_owner.values() for o in lst if o.type == family]
            before = [o for o in cands if o.id < op.id]
            pool = before or cands
            if pool:
                found = min(pool, key=lambda o: abs(op.id - o.id))
        if found is not None:
            self._claimed[id(op)] = found
        return found

    def _dispatch(self, op: Op):
        if op.block.get("IsSuppressed") == "true":
            return
        handler = getattr(self, "op_" + op.type.replace(" ", "_"), None)
        if handler is None:
            self.unsupported[op.type] = self.unsupported.get(op.type, 0) + 1
            self.warn(op.part_id, op.type, "operation not supported; skipped")
            return
        try:
            handler(op)
        except Exception as exc:  # noqa: BLE001 - keep replaying, report the failure
            self.warn(op.part_id, op.type, f"failed: {exc}", "error")

    # ------------------------------------------------------------------ primitives
    def _set(self, op: Op, shape):
        self.shapes[op.part_id] = self.to_global(shape, op.cs)

    def op_CreateRelativeCoordinateSystem(self, op):
        self._make_cs(op)

    def op_Point(self, op):
        pass

    def op_Box(self, op):
        q = op.params
        p0 = np.array([self.p(q, "XPosition"), self.p(q, "YPosition"), self.p(q, "ZPosition")])
        d = np.array([self.p(q, "XSize"), self.p(q, "YSize"), self.p(q, "ZSize")])
        lo, hi = np.minimum(p0, p0 + d), np.maximum(p0, p0 + d)
        self._set(op, BRepPrimAPI_MakeBox(_pnt(lo), _pnt(hi)).Shape())

    def op_Cylinder(self, op):
        q = op.params
        c = [self.p(q, "XCenter"), self.p(q, "YCenter"), self.p(q, "ZCenter")]
        r, h = self.p(q, "Radius"), self.p(q, "Height")
        axis = np.array(AXIS_VEC[q.get("WhichAxis", "Z")], float)
        if h < 0:
            axis, h = -axis, -h
        n = int(self.p(q, "NumSides", 0))
        if n >= 3:
            prof = self._ngon(c, r, axis, n)
            shape = BRepPrimAPI_MakePrism(prof, gp_Vec(*(axis * h))).Shape()
        else:
            shape = BRepPrimAPI_MakeCylinder(gp_Ax2(_pnt(c), gp_Dir(*axis)), r, h).Shape()
        self._set(op, shape)

    def op_Cone(self, op):
        q = op.params
        c = [self.p(q, "XCenter"), self.p(q, "YCenter"), self.p(q, "ZCenter")]
        r1, r2, h = self.p(q, "BottomRadius"), self.p(q, "TopRadius"), self.p(q, "Height")
        axis = np.array(AXIS_VEC[q.get("WhichAxis", "Z")], float)
        if h < 0:
            axis, h = -axis, -h
        self._set(op, BRepPrimAPI_MakeCone(gp_Ax2(_pnt(c), gp_Dir(*axis)), r1, r2, h).Shape())

    def op_Sphere(self, op):
        q = op.params
        c = [self.p(q, "XCenter"), self.p(q, "YCenter"), self.p(q, "ZCenter")]
        self._set(op, BRepPrimAPI_MakeSphere(_pnt(c), self.p(q, "Radius")).Shape())

    def op_Rectangle(self, op):
        q = op.params
        p0 = np.array([self.p(q, "XStart"), self.p(q, "YStart"), self.p(q, "ZStart")])
        a1, a2 = PLANE_AXES[q.get("WhichAxis", "Z")]
        u = np.array(AXIS_VEC[a1], float) * self.p(q, "Width")
        v = np.array(AXIS_VEC[a2], float) * self.p(q, "Height")
        poly = BRepBuilderAPI_MakePolygon(_pnt(p0), _pnt(p0 + u), _pnt(p0 + u + v), _pnt(p0 + v), True)
        self._set(op, poly.Wire())

    def op_Circle(self, op):
        q = op.params
        c = [self.p(q, "XCenter"), self.p(q, "YCenter"), self.p(q, "ZCenter")]
        r = self.p(q, "Radius")
        axis = np.array(AXIS_VEC[q.get("WhichAxis", "Z")], float)
        n = int(self.p(q, "NumSegments", 0))
        if n >= 3:
            wire = self._ngon(c, r, axis, n, as_face=False)
        else:
            circ = GC_MakeCircle(gp_Ax2(_pnt(c), gp_Dir(*axis)), r).Value()
            wire = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire()
        self._set(op, wire)

    def op_Ellipse(self, op):
        q = op.params
        from OCP.GC import GC_MakeEllipse
        c = [self.p(q, "XCenter"), self.p(q, "YCenter"), self.p(q, "ZCenter")]
        r1 = self.p(q, "MajRadius")
        r2 = r1 * self.p(q, "Ratio")
        axis = AXIS_VEC[q.get("WhichAxis", "Z")]
        a1, _ = PLANE_AXES[q.get("WhichAxis", "Z")]
        big, small = max(r1, r2), min(r1, r2)
        xdir = AXIS_VEC[a1] if r1 >= r2 else AXIS_VEC[PLANE_AXES[q.get("WhichAxis", "Z")][1]]
        el = GC_MakeEllipse(gp_Ax2(_pnt(c), gp_Dir(*axis), gp_Dir(*xdir)), big, small).Value()
        self._set(op, BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(el).Edge()).Wire())

    def _ngon(self, c, r, axis, n, as_face=True):
        axis = np.asarray(axis, float)
        ref = np.array([1.0, 0, 0]) if abs(axis[0]) < 0.9 else np.array([0, 1.0, 0])
        e1 = np.cross(axis, ref)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(axis, e1)
        poly = BRepBuilderAPI_MakePolygon()
        for k in range(n):
            t = 2 * math.pi * k / n
            poly.Add(_pnt(np.asarray(c) + r * (math.cos(t) * e1 + math.sin(t) * e2)))
        poly.Close()
        w = poly.Wire()
        return BRepBuilderAPI_MakeFace(w, True).Face() if as_face else w

    def op_Polyline(self, op):
        q = op.params
        pts_b = q.child("PolylinePoints")
        pts = [np.array([self.p(b, "X"), self.p(b, "Y"), self.p(b, "Z")]) for b in pts_b.children_named("PLPoint")]
        segs = q.child("PolylineSegments")
        seg_blocks = segs.children_named("PLSegment") if segs is not None else []
        if not seg_blocks:
            seg_blocks = [None] * (len(pts) - 1)
        closed = q.get("IsPolylineClosed") == "true"
        wire_mk = BRepBuilderAPI_MakeWire()
        last_end = None
        for i, sb in enumerate(seg_blocks):
            kind = sb.get("SegmentType", "Line") if sb is not None else "Line"
            start = int(sb.get("StartIndex", i)) if sb is not None else i
            npts = int(sb.get("NoOfPoints", 2)) if sb is not None else 2
            seg_pts = pts[start:start + npts]
            if last_end is not None:
                seg_pts[0] = last_end
            edge = None
            if kind == "Line":
                if np.linalg.norm(seg_pts[-1] - seg_pts[0]) > 1e-12:
                    edge = BRepBuilderAPI_MakeEdge(_pnt(seg_pts[0]), _pnt(seg_pts[-1])).Edge()
            elif kind == "Arc" and len(seg_pts) >= 3:
                edge = BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(*map(_pnt, seg_pts[:3])).Value()).Edge()
            elif kind == "AngularArc":
                center = np.array([self.p(sb, "ArcCenterX"), self.p(sb, "ArcCenterY"), self.p(sb, "ArcCenterZ")])
                ang = math.radians(self.p(sb, "ArcAngle"))
                plane = sb.get("ArcPlane", "XY")
                normal = {"XY": (0, 0, 1), "YZ": (1, 0, 0), "ZX": (0, 1, 0), "XZ": (0, -1, 0)}.get(plane, (0, 0, 1))
                p0 = seg_pts[0]
                rot = gp_Trsf()
                rot.SetRotation(gp_Ax1(_pnt(center), gp_Dir(*normal)), ang)
                pm = _pnt(p0).Transformed(_rot(center, normal, ang / 2))
                pe = _pnt(p0).Transformed(rot)
                edge = BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(_pnt(p0), pm, pe).Value()).Edge()
                seg_pts = [p0, np.array([pe.X(), pe.Y(), pe.Z()])]
            elif kind == "Spline" and len(seg_pts) >= 2:
                arr = TColgp_HArray1OfPnt(1, len(seg_pts))
                for k, sp in enumerate(seg_pts):
                    arr.SetValue(k + 1, _pnt(sp))
                interp = GeomAPI_Interpolate(arr, False, 1e-9)
                interp.Perform()
                edge = BRepBuilderAPI_MakeEdge(interp.Curve()).Edge()
            else:
                self.warn(op.part_id, "Polyline", f"segment type '{kind}' approximated by a line")
                edge = BRepBuilderAPI_MakeEdge(_pnt(seg_pts[0]), _pnt(seg_pts[-1])).Edge()
            if edge is not None:
                wire_mk.Add(edge)
            last_end = seg_pts[-1]
        if closed and last_end is not None and np.linalg.norm(last_end - pts[0]) > 1e-12:
            wire_mk.Add(BRepBuilderAPI_MakeEdge(_pnt(last_end), _pnt(pts[0])).Edge())
        wire = wire_mk.Wire()
        xs = q.child("PolylineXSection")
        xtype = xs.get("XSectionType", "None") if xs is not None else "None"
        if xtype != "None":
            wire = self._sweep_xsection(op, wire, pts, xs, xtype)
        self._set(op, wire)

    def _sweep_xsection(self, op, path, pts, xs, xtype):
        w = self.p(xs, "XSectionWidth", 0.0)
        h = self.p(xs, "XSectionHeight", 0.0)
        if abs(w) < 1e-15 and (xtype == "Line" or abs(h) < 1e-15):
            return path  # zero-size cross-section: AEDT keeps the plain polyline
        if xtype in ("Rectangle", "Isosceles Trapezoid") and abs(h) < 1e-12:
            xtype = "Line"  # zero-height trace is a flat ribbon (a degenerate solid crashes the kernel)
        p0 = pts[0]
        t0 = pts[1] - pts[0] if len(pts) > 1 else np.array([1.0, 0, 0])
        t0 = t0 / np.linalg.norm(t0)
        # plane normal of the path ('Auto' orientation): Z unless the path is vertical
        normal = _path_normal(pts)
        wdir = np.cross(normal, t0)
        wdir /= np.linalg.norm(wdir)
        if xtype == "Line":
            prof = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(_pnt(p0 - wdir * w / 2), _pnt(p0 + wdir * w / 2)).Edge()).Wire()
        elif xtype in ("Rectangle", "Isosceles Trapezoid"):
            wt = self.p(xs, "XSectionTopWidth", 0.0) if xtype != "Rectangle" else w
            wt = wt or w
            a, b = p0 - normal * h / 2, p0 + normal * h / 2
            prof = BRepBuilderAPI_MakePolygon(_pnt(a - wdir * w / 2), _pnt(a + wdir * w / 2),
                                              _pnt(b + wdir * wt / 2), _pnt(b - wdir * wt / 2), True).Wire()
        elif xtype == "Circle":
            circ = GC_MakeCircle(gp_Ax2(_pnt(p0), gp_Dir(*t0)), w / 2).Value()
            prof = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire()
        else:
            self.warn(op.part_id, "Polyline", f"cross-section '{xtype}' not supported; kept centre line")
            return path
        if len(pts) == 2 and count(path, TopAbs_EDGE) == 1:
            sweep = BRepPrimAPI_MakePrism(prof, gp_Vec(*(pts[1] - pts[0])))
            shape = sweep.Shape()
            if xtype != "Line":
                shape = _cap(shape, prof, pts[1] - pts[0])
            return shape
        pipe = BRepOffsetAPI_MakePipeShell(path)
        pipe.SetMode(gp_Dir(*normal))
        pipe.SetTransitionMode(BRepBuilderAPI_RightCorner)
        pipe.Add(prof, False, False)
        pipe.Build()
        if not pipe.IsDone():
            self.warn(op.part_id, "Polyline", f"could not sweep the '{xtype}' cross-section; kept centre line")
            return path
        if xtype != "Line":
            pipe.MakeSolid()
        return pipe.Shape()

    def op_Region(self, op):
        self.region_params = op.params
        self.region_part = op.part_id
        self.region_cs = op.cs

    def op_NativeBody(self, op):
        q = op.params
        src = q.get("SourceFile", "")
        shape = self._load_step(src)
        if shape is None:
            self.warn(op.part_id, "NativeBody",
                      f"imported body from '{os.path.basename(src) or '?'}' is stored in Parasolid format, "
                      "which cannot be read. Provide the original STEP file next to the .aedt to import it. "
                      "Faces/edges taken from it are rebuilt from their stored outlines where possible.")
            return
        self.shapes[op.part_id] = shape

    def _load_step(self, src):
        if not src:
            return None
        base = os.path.basename(src.replace("\\\\", "\\").replace("\\", "/"))
        candidates = [src] + [os.path.join(d, base) for d in self.step_dirs]
        if self.aedt_path:
            candidates.append(os.path.join(os.path.dirname(self.aedt_path), base))
        for c in candidates:
            if c and os.path.exists(c) and os.path.splitext(c)[1].lower() in (".step", ".stp"):
                if c not in self._step_cache:
                    r = STEPControl_Reader()
                    if r.ReadFile(c) != 1:
                        continue
                    r.TransferRoots()
                    shape = r.OneShape()
                    scale = {"mm": 1.0}.get(self.units, None)
                    if scale is None:
                        self.issues.append(Issue("warn", "", "NativeBody", "STEP assumed to be in model units"))
                    self._step_cache[c] = shape
                return self._step_cache[c]
        return None

    # ------------------------------------------------------------------ local operations
    def op_CoverLines(self, op):
        shape = self.shapes.get(op.part_id)
        if is_null(shape):
            raise ValueError("nothing to cover")
        faces = []
        for w in _iter(shape, TopAbs_WIRE):
            mk = BRepBuilderAPI_MakeFace(TopoDS.Wire(w), True)
            if mk.IsDone():
                faces.append(mk.Face())
        if not faces:
            raise ValueError("lines are not closed/planar")
        # nested wires (holes) are handled by subtracting inner faces from outer ones
        faces.sort(key=face_area, reverse=True)
        result = faces[0]
        for f in faces[1:]:
            if dist(f, result) < 1e-9 and face_area(f) < face_area(faces[0]):
                result = BRepAlgoAPI_Cut(result, f).Shape()
            else:
                result = compound([result, f])
        self.shapes[op.part_id] = result

    def op_CoverSurfaces(self, op):
        self.warn(op.part_id, "CoverSurfaces", "treated as solid")

    def op_Move(self, op):
        q = op.params
        v = [self.p(q, "TranslateVectorX"), self.p(q, "TranslateVectorY"), self.p(q, "TranslateVectorZ")]
        v = self.cs_vec(v, op.cs)
        t = gp_Trsf()
        t.SetTranslation(gp_Vec(*v))
        self._apply(op, t)

    def op_Rotate(self, op):
        q = op.params
        axis = self.cs_vec(AXIS_VEC[q.get("RotateAxis", "Z")], op.cs)
        origin = self.cs_point((0, 0, 0), op.cs)
        t = _rot(origin, axis, math.radians(self.p(q, "RotateAngle")))
        self._apply(op, t)

    def op_Mirror(self, op):
        q = op.params
        base = self.cs_point([self.p(q, "MirrorBase" + c) for c in "XYZ"], op.cs)
        n = self.cs_vec([self.p(q, "MirrorNormal" + c) for c in "XYZ"], op.cs)
        t = gp_Trsf()
        t.SetMirror(gp_Ax2(_pnt(base), gp_Dir(*n)))
        self._apply(op, t)

    def op_Scale(self, op):
        q = op.params
        f = [self.p(q, k, 1.0) for k in ("ScaleX", "ScaleY", "ScaleZ")]
        if max(f) - min(f) > 1e-12:
            raise ValueError("non-uniform scaling not supported")
        t = gp_Trsf()
        t.SetScale(_pnt(self.cs_point((0, 0, 0), op.cs)), f[0])
        self._apply(op, t)

    def _apply(self, op, t):
        target = int(op.params.get("TargetID", op.part_id)) if op.params is not None else op.part_id
        for pid in (op.part_id, target):
            if not is_null(self.shapes.get(pid)):
                self.shapes[pid] = transformed(self.shapes[pid], t)
                acc = gp_Trsf()
                acc.Multiply(t)
                acc.Multiply(self.part_trsf.get(pid, gp_Trsf()))
                self.part_trsf[pid] = acc
                return
        raise ValueError("target has no geometry")

    def _boolean(self, op, kind):
        q = op.params
        blank = int(q.get("BlankPart", op.part_id))
        tools = _idlist(q, "ToolParts")
        keep = q.get("KeepOriginals") == "true"
        shape = self.shapes.get(blank)
        tool_shapes = [self.shapes.get(t) for t in tools]
        missing = [t for t, s in zip(tools, tool_shapes) if is_null(s)]
        if missing:
            self.warn(blank, op.type, f"{len(missing)} tool object(s) have no geometry; ignored")
        tool_shapes = [s for s in tool_shapes if not is_null(s)]
        if is_null(shape):
            if kind == "fuse" and tool_shapes:
                self.warn(blank, op.type, "blank object has no geometry; result uses tool objects only")
                shape = tool_shapes.pop(0)
            else:
                raise ValueError("blank object has no geometry")
        algo = {"fuse": BRepAlgoAPI_Fuse, "cut": BRepAlgoAPI_Cut, "common": BRepAlgoAPI_Common}[kind]
        for t in tool_shapes:
            mk = algo(shape, t)
            result = mk.Shape() if mk.IsDone() else None
            if is_null(result) or (kind != "common" and count(result, TopAbs_VERTEX) == 0):
                # The kernel could not combine these (typically a solid with a sheet).
                if kind == "fuse":
                    result = compound([shape, t])
                elif kind == "cut":
                    self.warn(blank, op.type, "kernel could not subtract one tool; it was ignored")
                    result = shape
                else:
                    self.warn(blank, op.type, "kernel could not intersect; kept the blank object")
                    result = shape
            shape = result
        self.shapes[blank] = shape
        if not keep:
            for t in tools:
                self.shapes.pop(t, None)

    def op_Unite(self, op):
        self._boolean(op, "fuse")

    def op_Substract(self, op):
        self._boolean(op, "cut")

    op_Subtract = op_Substract

    def op_Intersect(self, op):
        self._boolean(op, "common")

    def op_SweepAlongVector(self, op):
        q = op.params
        pid = int(q.get("ProfileID", op.part_id))
        v = self.cs_vec([self.p(q, "SweepVector" + c) for c in "XYZ"], op.cs)
        shape = self.shapes.get(pid)
        if is_null(shape):
            raise ValueError("profile has no geometry")
        if abs(self.p(q, "DraftAngle", 0.0)) > 1e-9:
            self.warn(pid, op.type, "draft angle ignored")
        self.shapes[pid] = BRepPrimAPI_MakePrism(shape, gp_Vec(*v)).Shape()

    def op_SweepAroundAxis(self, op):
        from OCP.BRepPrimAPI import BRepPrimAPI_MakeRevol
        q = op.params
        pid = int(q.get("ProfileID", op.part_id))
        axis = self.cs_vec(AXIS_VEC[q.get("SweepAxis", "Z")], op.cs)
        ang = math.radians(self.p(q, "SweepAngle"))
        shape = self.shapes.get(pid)
        self.shapes[pid] = BRepPrimAPI_MakeRevol(shape, gp_Ax1(_pnt(self.cs_point((0, 0, 0), op.cs)), gp_Dir(*axis)), ang).Shape()

    def op_SheetThicken(self, op):
        q = op.params
        pid = int(q.get("LocalOpPart", op.part_id))
        t = self.p(q, "Thickness")
        both = q.get("BothSides") == "true"
        shape = self.shapes.get(pid)
        if is_null(shape):
            raise ValueError("sheet has no geometry")
        sign = 1.0
        info = q.find("ThickenAdditionalInfo", "ShellThickenDirectionInfo")
        solids = []
        for f in faces_of(shape):
            n = face_normal(f)
            if n is None:
                raise ValueError("only planar sheets can be thickened")
            if info is not None and info.child("DirectionAtPoint") is not None:
                d = [self.p(info.child("DirectionAtPoint"), c, 0.0) for c in "XYZ"]
                if np.dot(d, n) < 0:
                    sign = -1.0
            if both:
                tr = gp_Trsf()
                tr.SetTranslation(gp_Vec(*(-n * t)))
                f = transformed(f, tr)
                solids.append(BRepPrimAPI_MakePrism(f, gp_Vec(*(n * 2 * t))).Shape())
            else:
                solids.append(BRepPrimAPI_MakePrism(f, gp_Vec(*(n * t * sign))).Shape())
        self.shapes[pid] = compound(solids)

    def op_Fillet(self, op):
        self._blend(op, fillet=True)

    def op_Chamfer(self, op):
        self._blend(op, fillet=False)

    def _blend(self, op, fillet):
        q = op.params
        pid = int(q.get("PartID", op.part_id))
        shape = self.shapes.get(pid)
        if is_null(shape) or count(shape, TopAbs_SOLID) == 0:
            raise ValueError("fillet/chamfer only supported on solids")
        edges = [self.resolve_edge(pid, e) for e in _idlist(q, "Edges")]
        edges = [e for e in edges if e is not None]
        if not edges:
            raise ValueError("could not locate the edges")
        if fillet:
            mk = BRepFilletAPI_MakeFillet(shape)
            r = self.p(q, "Radius")
            for e in edges:
                mk.Add(r, e)
        else:
            mk = BRepFilletAPI_MakeChamfer(shape)
            d = self.p(q, "LeftDistance")
            for e in edges:
                mk.Add(d, e)
        mk.Build()
        if not mk.IsDone():
            raise ValueError("kernel could not build the blend; sharp corner kept")
        self.shapes[pid] = mk.Shape()

    # ------------------------------------------------------------------ To/From pairs
    def op_BodyFromFace(self, op):
        q = op.params
        pid = int(q.get("LocalOpPart", op.part_id))
        ids = _idlist(q, "FacesToDetach", "Faces")
        faces = [self.resolve_face(pid, f) for f in ids]
        if any(f is None for f in faces):
            self.warn(pid, op.type, "some faces could not be located")
        self.pending[id(op)] = [BRepBuilderAPI_Copy(f).Shape() if f is not None else None for f in faces]

    def op_BodyFromFaceFrom(self, op):
        self._take(op, "BodyFromFace", "WhichFace")

    def op_BodyFromEdge(self, op):
        q = op.params
        pid = int(q.get("LocalOpPart", op.part_id))
        edges = [self.resolve_edge(pid, e) for e in _idlist(q, "Edges")]
        if any(e is None for e in edges):
            self.warn(pid, op.type, "some edges could not be located")
        self.pending[id(op)] = [BRepBuilderAPI_MakeWire(e).Wire() if e is not None else None for e in edges]

    def op_BodyFromEdgeFrom(self, op):
        self._take(op, "BodyFromEdge", "WhichEdge")

    def op_SweepFaceAlongNormal(self, op):
        q = op.params
        pid = int(q.get("LocalOpPart", op.part_id))
        length = self.p(q, "LengthOfSweep")
        out = []
        for fid in _idlist(q, "FacesToDetach"):
            f = self.resolve_face(pid, fid)
            n = face_normal(f) if f is not None else None
            if f is None or n is None:
                self.warn(pid, op.type, f"face {fid} could not be located")
                out.append(None)
                continue
            out.append(BRepPrimAPI_MakePrism(f, gp_Vec(*(n * length))).Shape())
        self.pending[id(op)] = out

    def op_SweepFaceAlongNormalFrom(self, op):
        self._take(op, "SweepFaceAlongNormal", "WhichFace")

    def op_CloneTo(self, op):
        # snapshot the source as it is at this point of its history
        shape = self.shapes.get(op.owner)
        self.pending[id(op)] = [BRepBuilderAPI_Copy(shape).Shape() if not is_null(shape) else None]

    def op_CloneFrom(self, op):
        to_op = self._claimed.get(id(op))
        if to_op is not None and self.pending.get(id(to_op)):
            snap = self.pending[id(to_op)][0]
            if snap is None:
                raise ValueError("clone source has no geometry")
            self.shapes[op.part_id] = BRepBuilderAPI_Copy(snap).Shape()
            return
        src = int(op.params.get("SourceID"))
        shape = self.shapes.get(src)
        if is_null(shape):
            raise ValueError("clone source has no geometry")
        self.shapes[op.part_id] = BRepBuilderAPI_Copy(shape).Shape()

    def op_DuplicateMirror(self, op):
        q = op.params
        base = self.cs_point([self.p(q, "DuplicateMirrorBase" + c) for c in "XYZ"], op.cs)
        n = self.cs_vec([self.p(q, "DuplicateMirrorNormal" + c) for c in "XYZ"], op.cs)
        t = gp_Trsf()
        t.SetMirror(gp_Ax2(_pnt(base), gp_Dir(*n)))
        shape = self.shapes.get(op.part_id)
        self.pending[id(op)] = [transformed(shape, t) if not is_null(shape) else None]

    def op_DuplicateBodyMirror(self, op):
        self._take(op, "DuplicateMirror", "WhichClone")

    def op_DuplicateAlongLine(self, op):
        q = op.params
        v = self.cs_vec([self.p(q, "XComponent"), self.p(q, "YComponent"), self.p(q, "ZComponent")], op.cs)
        n = int(self.p(q, "NumClones", 2))
        shape = self.shapes.get(op.part_id)
        out = []
        for k in range(1, n):
            t = gp_Trsf()
            t.SetTranslation(gp_Vec(*(np.asarray(v) * k)))
            out.append(transformed(shape, t))
        self.pending[id(op)] = out

    def op_DuplicateBodyAlongLine(self, op):
        self._take(op, "DuplicateAlongLine", "WhichClone")

    def op_DuplicateAroundAxis(self, op):
        q = op.params
        axis = self.cs_vec(AXIS_VEC[q.get("WhichAxis", "Z")], op.cs)
        ang = math.radians(self.p(q, "AngleStr"))
        n = int(self.p(q, "NumClones", 2))
        shape = self.shapes.get(op.part_id)
        origin = self.cs_point((0, 0, 0), op.cs)
        self.pending[id(op)] = [transformed(shape, _rot(origin, axis, ang * k)) for k in range(1, n)]

    def op_DuplicateBodyAroundAxis(self, op):
        self._take(op, "DuplicateAroundAxis", "WhichClone")

    def _take(self, op, family, which_key):
        to_op = self._claimed.get(id(op))
        items = (self.pending.get(id(to_op)) if to_op is not None else None) or []
        k = int(op.params.get(which_key, "0")) if op.params is not None else 0
        if k >= len(items) or items[k] is None:
            raise ValueError(f"source geometry for {family} unavailable")
        self.shapes[op.part_id] = items[k]

    # ------------------------------------------------------------------ region
    def _region_box(self):
        q = self.region_params
        model = [s for pid, s in self.shapes.items()
                 if pid != self.region_part and pid in self.parts and self.parts[pid].toplevel
                 and self.parts[pid].model and not is_null(s)]
        if not model:
            return None
        lo = np.full(3, np.inf)
        hi = np.full(3, -np.inf)
        for s in model:
            a, b = bounds(s)
            lo, hi = np.minimum(lo, a), np.maximum(hi, b)
        size = hi - lo
        out_lo, out_hi = lo.copy(), hi.copy()
        for i, ax in enumerate("XYZ"):
            for sign, arr in (("+", out_hi), ("-", out_lo)):
                ptype = q.get(f"{sign}{ax}PaddingType", "Percentage Offset")
                val = q.get(f"{sign}{ax}Padding", "0")
                if ptype == "Percentage Offset":
                    off = float(self.num(val)) / 100 * size[i]
                elif ptype == "Absolute Position":
                    arr[i] = self.num(val)
                    continue
                else:
                    off = self.num(val)
                arr[i] += off if sign == "+" else -off
        return {"lo": out_lo, "hi": out_hi}


# --------------------------------------------------------------------------- utilities

def _idlist(q: Block, *names) -> list[int]:
    """Integer id list stored either as an array entry (``Edges[1: 715]``) or a property."""
    for n in names:
        c = q.call(n)
        if c is not None:
            return ints(c)
        if n in q.props:
            return ints(q.props[n])
    return []


def _rot(origin, axis, angle) -> gp_Trsf:
    t = gp_Trsf()
    t.SetRotation(gp_Ax1(_pnt(origin), gp_Dir(*map(float, axis))), angle)
    return t


def _path_normal(pts) -> np.ndarray:
    pts = np.asarray(pts, float)
    if len(pts) >= 3:
        c = pts.mean(axis=0)
        _, s, vt = np.linalg.svd(pts - c)
        if s[1] > 1e-9 * max(s[0], 1e-12):
            n = vt[2]
            return n if n[np.argmax(np.abs(n))] > 0 else -n
    t = pts[-1] - pts[0]
    t /= np.linalg.norm(t)
    return np.array([0, 0, 1.0]) if abs(t[2]) < 0.9 else np.array([1.0, 0, 0])


def _cap(shell, prof, vec):
    """Close an extruded tube (from a closed profile wire) into a solid."""
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid
    f0 = BRepBuilderAPI_MakeFace(prof, True).Face()
    t = gp_Trsf()
    t.SetTranslation(gp_Vec(*vec))
    f1 = transformed(f0, t)
    sew = BRepBuilderAPI_Sewing(1e-7)
    for s in (shell, f0, f1):
        sew.Add(s)
    sew.Perform()
    sewn = sew.SewedShape()
    shells = [TopoDS.Shell(s) for s in _iter(sewn, TopAbs_SHELL)]
    if not shells:
        return BRepPrimAPI_MakePrism(f0, gp_Vec(*vec)).Shape()
    mk = BRepBuilderAPI_MakeSolid(shells[0])
    return mk.Solid() if mk.IsDone() else BRepPrimAPI_MakePrism(f0, gp_Vec(*vec)).Shape()


def _fit_circle(pts2d):
    """Least-squares circle through 2D points -> (centre, radius, max residual)."""
    x, y = pts2d[:, 0], pts2d[:, 1]
    A = np.column_stack([x, y, np.ones_like(x)])
    sol, *_ = np.linalg.lstsq(A, x ** 2 + y ** 2, rcond=None)
    cx, cy = sol[0] / 2, sol[1] / 2
    r = math.sqrt(max(sol[2] + cx ** 2 + cy ** 2, 0.0))
    res = np.abs(np.hypot(x - cx, y - cy) - r).max()
    return np.array([cx, cy]), r, res


def _outline_wire(verts, origin, e1, e2, with_arcs=True):
    """Closed wire through outline vertices.

    AEDT stores only the vertices of a face, so an arc edge (two end vertices,
    no interior points) would come out as a straight chord. Where an unusually
    long chord is flanked by vertices lying on a common circle, it is replaced
    by the arc of that circle.
    """
    v = np.asarray(verts, float)
    n = len(v)
    uv = np.column_stack([(v - origin) @ e1, (v - origin) @ e2])
    seg = np.linalg.norm(np.roll(v, -1, axis=0) - v, axis=1)
    med = float(np.median(seg))
    mk = BRepBuilderAPI_MakeWire()
    for i in range(n):
        a, b = v[i], v[(i + 1) % n]
        if seg[i] < 1e-12:
            continue
        edge = None
        if with_arcs and n >= 10 and seg[i] > 4 * med:
            idx = [(i - k) % n for k in range(0, 6)] + [(i + 1 + k) % n for k in range(0, 6)]
            centre, r, res = _fit_circle(uv[idx])
            if r > 0 and res < 2e-3 * r and r < 1e3 * seg[i]:
                ua, ub = uv[i], uv[(i + 1) % n]
                ta = math.atan2(*(ua - centre)[::-1])
                tb = math.atan2(*(ub - centre)[::-1])
                # take the arc that does not run through the neighbouring vertices
                span = (tb - ta) % (2 * math.pi)
                nb = uv[(i - 2) % n] - centre
                tn = (math.atan2(nb[1], nb[0]) - ta) % (2 * math.pi)
                tm = ta + span / 2 if not (0 < tn < span) else ta + span / 2 + math.pi
                m2 = centre + r * np.array([math.cos(tm), math.sin(tm)])
                m3 = origin + m2[0] * e1 + m2[1] * e2
                try:
                    edge = BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(_pnt(a), _pnt(m3), _pnt(b)).Value()).Edge()
                except Exception:  # noqa: BLE001
                    edge = None
        if edge is None:
            edge = BRepBuilderAPI_MakeEdge(_pnt(a), _pnt(b)).Edge()
        mk.Add(edge)
    return mk.Wire() if mk.IsDone() else None


def _order_loop(verts, e1, e2):
    """Order unsorted boundary vertices of a planar face into a loop (greedy nearest neighbour)."""
    pts = np.asarray(verts)
    n = len(pts)
    remaining = list(range(1, n))
    order = [0]
    while remaining:
        last = pts[order[-1]]
        k = min(remaining, key=lambda j: np.sum((pts[j] - last) ** 2))
        order.append(k)
        remaining.remove(k)
    return pts[order]


def _part_info(pid, attrs: Block | None, toplevel: bool) -> PartInfo:
    if attrs is None:
        return PartInfo(pid, f"Part{pid}", toplevel=toplevel)
    color = (0.6, 0.6, 0.6)
    c = attrs.get("Color")
    if c:
        try:
            r, g, b = (int(v) for v in c.strip("()").split())
            color = (r / 255, g / 255, b / 255)
        except ValueError:
            pass
    flags = attrs.get("Flags", "") or ""
    return PartInfo(
        pid, attrs.get("Name", f"Part{pid}"),
        material=(attrs.get("MaterialValue", "") or "").strip('"'),
        solve_inside=attrs.get("SolveInside", "true") == "true",
        color=color, transparency=float(attrs.get("Transparency", "0") or 0),
        model="NonModel" not in flags, toplevel=toplevel,
    )


# --------------------------------------------------------------------------- tessellation

def tessellate(shape, deflection: float) -> tuple[np.ndarray, np.ndarray]:
    """Triangulate a shape; returns (points Nx3, triangles Mx3) with coincident nodes welded."""
    BRepMesh_IncrementalMesh(shape, deflection, False, 0.3, True)
    pts, tris, off = [], [], 0
    for f in faces_of(shape):
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(f, loc)
        if tri is None:
            continue
        trsf = loc.Transformation()
        nodes = np.array([[p.X(), p.Y(), p.Z()] for p in
                          (tri.Node(i).Transformed(trsf) for i in range(1, tri.NbNodes() + 1))])
        t = np.array([tri.Triangle(i).Get() for i in range(1, tri.NbTriangles() + 1)]) - 1
        if f.Orientation() == TopAbs_REVERSED:
            t = t[:, ::-1]
        pts.append(nodes)
        tris.append(t + off)
        off += len(nodes)
    if not pts:
        return np.zeros((0, 3)), np.zeros((0, 3), int)
    pts = np.vstack(pts)
    tris = np.vstack(tris)
    # weld duplicated boundary nodes so closed solids become watertight
    scale = max(np.ptp(pts, axis=0).max(), 1e-12)
    key = np.round(pts / (scale * 1e-9)).astype(np.int64)
    _, idx, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    pts = pts[idx]
    tris = inv.reshape(-1)[tris]
    good = (tris[:, 0] != tris[:, 1]) & (tris[:, 1] != tris[:, 2]) & (tris[:, 0] != tris[:, 2])
    return pts, tris[good]


def shape_kind(shape) -> str:
    if count(shape, TopAbs_SOLID):
        return "solid"
    if count(shape, TopAbs_FACE):
        return "sheet"
    if count(shape, TopAbs_EDGE):
        return "wire"
    return "empty"


def sheet_to_slab(shape, thickness: float):
    """Give zero-thickness sheets a tiny thickness so they become closed polyhedra for the FDTD mesher.

    Adjacent coplanar faces are merged first; thickening them one by one would leave shared side walls
    and produce a non-manifold (not watertight) polyhedron.
    """
    from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
    try:
        fused = shape
        fs = faces_of(shape)
        if len(fs) > 1:
            fused = fs[0]
            for f in fs[1:]:
                fused = BRepAlgoAPI_Fuse(fused, f).Shape()
        unify = ShapeUpgrade_UnifySameDomain(fused, True, True, True)
        unify.Build()
        shape = unify.Shape()
    except Exception:  # noqa: BLE001 - fall back to the faces as they are
        pass
    out = []
    for f in faces_of(shape):
        n = face_normal(f)
        if n is None:
            continue
        th = thickness
        if callable(thickness):
            # per-face thickness: (normal axis, plane coordinate) for axis-aligned faces, else (None, None)
            ax = int(np.argmax(np.abs(n)))
            if abs(abs(n[ax]) - 1) < 1e-9:
                lo, hi = bounds(f)
                th = thickness(ax, 0.5 * (lo[ax] + hi[ax]))
            else:
                th = thickness(None, None)
        t = gp_Trsf()
        t.SetTranslation(gp_Vec(*(-n * th / 2)))
        out.append(BRepPrimAPI_MakePrism(transformed(f, t), gp_Vec(*(n * th))).Shape())
    return compound(out)


def solids_of(shape):
    return [TopoDS.Solid(s) for s in _iter(shape, TopAbs_SOLID)]


def write_brep(shape, path: str):
    from OCP.BRepTools import BRepTools
    BRepTools.Write_s(shape, path)


def read_brep(path: str):
    from OCP.BRepTools import BRepTools
    shape = TopoDS_Shape()
    BRepTools.Read_s(shape, path, BRep_Builder())
    return shape
