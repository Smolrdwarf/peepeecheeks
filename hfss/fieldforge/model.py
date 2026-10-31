"""Project data model: geometry, materials, ports, boundaries, solution setup.

Every numeric geometry field is stored as an expression string so that
design variables (e.g. ``L = 28.5``) can be referenced and swept.
"""
from __future__ import annotations

import ast
import json
import math
import operator
from dataclasses import asdict, dataclass, field

FORMAT_VERSION = 1

UNITS = {"m": 1.0, "mm": 1e-3, "um": 1e-6, "mil": 25.4e-6, "in": 25.4e-3}
AXES = ("X", "Y", "Z")


# --------------------------------------------------------------------------- expressions

_BINOPS = {
    ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
    ast.Div: operator.truediv, ast.Pow: operator.pow, ast.Mod: operator.mod,
}
_UNARY = {ast.UAdd: operator.pos, ast.USub: operator.neg}
_FUNCS = {
    "sin": math.sin, "cos": math.cos, "tan": math.tan, "sqrt": math.sqrt,
    "abs": abs, "min": min, "max": max, "exp": math.exp, "log": math.log,
    "log10": math.log10, "atan": math.atan, "asin": math.asin, "acos": math.acos,
}
_CONSTS = {"pi": math.pi, "c0": 299792458.0, "e": math.e}


class ExpressionError(ValueError):
    pass


def evaluate(expr, variables: dict[str, float] | None = None) -> float:
    """Safely evaluate an arithmetic expression that may reference variables."""
    if isinstance(expr, (int, float)):
        return float(expr)
    variables = variables or {}
    try:
        tree = ast.parse(str(expr).strip(), mode="eval")
    except SyntaxError as exc:
        raise ExpressionError(f"Invalid expression '{expr}'") from exc

    def ev(node):
        if isinstance(node, ast.Expression):
            return ev(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.BinOp) and type(node.op) in _BINOPS:
            return _BINOPS[type(node.op)](ev(node.left), ev(node.right))
        if isinstance(node, ast.UnaryOp) and type(node.op) in _UNARY:
            return _UNARY[type(node.op)](ev(node.operand))
        if isinstance(node, ast.Name):
            if node.id in variables:
                return float(variables[node.id])
            if node.id in _CONSTS:
                return _CONSTS[node.id]
            raise ExpressionError(f"Unknown variable '{node.id}' in '{expr}'")
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in _FUNCS:
            return float(_FUNCS[node.func.id](*[ev(a) for a in node.args]))
        raise ExpressionError(f"Unsupported syntax in '{expr}'")

    return ev(tree)


# --------------------------------------------------------------------------- materials

@dataclass
class Material:
    name: str
    epsr: float = 1.0
    mur: float = 1.0
    tand: float = 0.0           # dielectric loss tangent
    conductor: bool = False     # modelled as PEC in the FDTD solver
    color: tuple = (0.6, 0.6, 0.6)


MATERIAL_LIBRARY = {m.name: m for m in [
    Material("vacuum", color=(0.75, 0.85, 1.0)),
    Material("air", 1.0006, color=(0.75, 0.85, 1.0)),
    Material("pec", conductor=True, color=(0.85, 0.65, 0.2)),
    Material("copper", conductor=True, color=(0.85, 0.5, 0.3)),
    Material("aluminum", conductor=True, color=(0.78, 0.78, 0.82)),
    Material("FR4_epoxy", 4.4, tand=0.02, color=(0.25, 0.6, 0.3)),
    Material("Rogers RO4003C", 3.55, tand=0.0027, color=(0.85, 0.8, 0.55)),
    Material("Rogers RT/duroid 5880", 2.2, tand=0.0009, color=(0.9, 0.9, 0.8)),
    Material("Teflon (PTFE)", 2.1, tand=0.0002, color=(0.95, 0.95, 0.95)),
    Material("Alumina_96pct", 9.4, tand=0.0006, color=(0.95, 0.92, 0.85)),
    Material("silicon", 11.9, color=(0.4, 0.4, 0.5)),
    Material("glass", 5.5, tand=0.005, color=(0.6, 0.85, 0.9)),
    Material("water_distilled", 81.0, tand=0.05, color=(0.3, 0.5, 0.9)),
]}


# --------------------------------------------------------------------------- geometry

# Parameter schemas per primitive kind: (key, label, default expression)
PRIMITIVE_PARAMS = {
    "box": [("x", "Position X", "0"), ("y", "Position Y", "0"), ("z", "Position Z", "0"),
            ("dx", "X Size", "10"), ("dy", "Y Size", "10"), ("dz", "Z Size", "10")],
    "cylinder": [("x", "Center X", "0"), ("y", "Center Y", "0"), ("z", "Center Z", "0"),
                 ("radius", "Radius", "5"), ("height", "Height", "10")],
    "sphere": [("x", "Center X", "0"), ("y", "Center Y", "0"), ("z", "Center Z", "0"),
               ("radius", "Radius", "5")],
    "rectangle": [("x", "Position X", "0"), ("y", "Position Y", "0"), ("z", "Position Z", "0"),
                  ("w", "Width", "10"), ("l", "Length", "10")],
}
# Kinds that carry an orientation axis (cylinder axis / rectangle normal)
AXIS_KINDS = {"cylinder", "rectangle"}


def _in_plane_axes(normal: int) -> tuple[int, int]:
    """Cyclic in-plane axes for a sheet with the given normal (Z -> X,Y; X -> Y,Z; Y -> Z,X)."""
    return (normal + 1) % 3, (normal + 2) % 3


@dataclass
class Primitive:
    name: str
    kind: str
    material: str = "vacuum"
    params: dict = field(default_factory=dict)
    axis: str = "Z"
    transparency: float = 0.0
    visible: bool = True

    def __post_init__(self):
        for key, _, default in PRIMITIVE_PARAMS[self.kind]:
            self.params.setdefault(key, default)

    def values(self, variables) -> dict[str, float]:
        return {k: evaluate(v, variables) for k, v in self.params.items()}

    def bounds(self, variables) -> tuple[list[float], list[float]]:
        """Axis-aligned bounding box (start, stop) in project units."""
        v = self.values(variables)
        c = [v["x"], v["y"], v["z"]]
        a = AXES.index(self.axis)
        if self.kind == "box":
            p2 = [c[0] + v["dx"], c[1] + v["dy"], c[2] + v["dz"]]
            return [min(c[i], p2[i]) for i in range(3)], [max(c[i], p2[i]) for i in range(3)]
        if self.kind == "sphere":
            r = v["radius"]
            return [ci - r for ci in c], [ci + r for ci in c]
        if self.kind == "cylinder":
            r = v["radius"]
            lo = [ci - r for ci in c]
            hi = [ci + r for ci in c]
            lo[a], hi[a] = sorted((c[a], c[a] + v["height"]))
            return lo, hi
        # rectangle
        u, w = _in_plane_axes(a)
        p2 = list(c)
        p2[u] += v["w"]
        p2[w] += v["l"]
        return [min(c[i], p2[i]) for i in range(3)], [max(c[i], p2[i]) for i in range(3)]


@dataclass
class LumpedPort:
    name: str
    # Start/stop corners of the port sheet; current flows along `axis`.
    params: dict = field(default_factory=lambda: {
        "x1": "0", "y1": "0", "z1": "0", "x2": "0", "y2": "0", "z2": "1"})
    axis: str = "Z"
    impedance: str = "50"
    excited: bool = True

    def corners(self, variables):
        v = {k: evaluate(e, variables) for k, e in self.params.items()}
        return [v["x1"], v["y1"], v["z1"]], [v["x2"], v["y2"], v["z2"]]


BOUNDARY_TYPES = ["PML_8", "MUR", "PEC", "PMC"]
FACE_NAMES = ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]


@dataclass
class Region:
    """Simulation domain around the model, with a boundary condition per face."""
    padding: str = "auto"                                   # 'auto' = lambda/4 at f_min
    boundaries: list = field(default_factory=lambda: ["PML_8"] * 6)


@dataclass
class Setup:
    name: str = "Setup1"
    f_start: str = "1"          # GHz
    f_stop: str = "3"           # GHz
    n_points: int = 401
    cells_per_wavelength: int = 20
    end_criteria_db: float = -40.0
    max_timesteps: int = 60000
    farfield: bool = True
    farfield_freq: str = "2"    # GHz
    farfield_step_deg: float = 3.0


@dataclass
class Project:
    name: str = "Project1"
    units: str = "mm"
    variables: dict = field(default_factory=dict)          # name -> expression
    objects: list = field(default_factory=list)
    ports: list = field(default_factory=list)
    materials: dict = field(default_factory=dict)          # user-defined, name -> Material
    region: Region = field(default_factory=Region)
    setup: Setup = field(default_factory=Setup)

    # ----------------------------------------------------------------- helpers
    def material(self, name: str) -> Material:
        if name in self.materials:
            return self.materials[name]
        return MATERIAL_LIBRARY.get(name, MATERIAL_LIBRARY["vacuum"])

    def all_materials(self) -> dict[str, Material]:
        return {**MATERIAL_LIBRARY, **self.materials}

    def resolved_variables(self) -> dict[str, float]:
        """Evaluate design variables in order; later variables may use earlier ones."""
        out: dict[str, float] = {}
        pending = dict(self.variables)
        for _ in range(len(pending) + 1):
            progressed = False
            for name, expr in list(pending.items()):
                try:
                    out[name] = evaluate(expr, out)
                except ExpressionError:
                    continue
                del pending[name]
                progressed = True
            if not pending or not progressed:
                break
        if pending:
            name, expr = next(iter(pending.items()))
            evaluate(expr, out)  # re-raise with a useful message
        return out

    def unique_name(self, base: str) -> str:
        taken = {o.name for o in self.objects} | {p.name for p in self.ports}
        i = 1
        while f"{base}{i}" in taken:
            i += 1
        return f"{base}{i}"

    def unit_scale(self) -> float:
        return UNITS[self.units]

    # ----------------------------------------------------------------- (de)serialisation
    def to_dict(self) -> dict:
        d = asdict(self)
        d["format_version"] = FORMAT_VERSION
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "Project":
        d = dict(d)
        d.pop("format_version", None)
        p = cls(
            name=d.get("name", "Project1"),
            units=d.get("units", "mm"),
            variables=dict(d.get("variables", {})),
            objects=[Primitive(**o) for o in d.get("objects", [])],
            ports=[LumpedPort(**o) for o in d.get("ports", [])],
            materials={k: Material(**{**m, "color": tuple(m.get("color", (0.6, 0.6, 0.6)))})
                       for k, m in d.get("materials", {}).items()},
            region=Region(**d.get("region", {})),
            setup=Setup(**d.get("setup", {})),
        )
        return p

    def save(self, path: str) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=2)

    @classmethod
    def load(cls, path: str) -> "Project":
        with open(path, encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))
