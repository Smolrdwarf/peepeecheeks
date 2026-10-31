"""Translate a Project into an openEMS FDTD simulation and post-process it.

Run as a subprocess so the GUI stays responsive and the solve can be aborted:

    python -m fieldforge.solver <project.json> <output_dir>

Results are written to ``<output_dir>/results.npz``.
"""
from __future__ import annotations

import os
import sys

import numpy as np

from .model import AXES, Project, evaluate

C0 = 299792458.0
EPS0 = 8.8541878128e-12
# Below this fraction of accepted power (1 - |S11|^2) radiated/accepted is dominated by numerical noise
MIN_ACCEPTED = 0.01


def openems_install_path() -> str:
    # Windows uses the release zip; macOS/Linux build from source into the upstream default prefix.
    default = ("tools", "openEMS") if os.name == "nt" else ("opt", "openEMS")
    return os.environ.get("OPENEMS_INSTALL_PATH", os.path.join(os.path.expanduser("~"), *default))


def _import_openems():
    path = openems_install_path()
    os.environ["OPENEMS_INSTALL_PATH"] = path
    if os.name == "nt" and os.path.isdir(path):
        os.add_dll_directory(path)
        os.environ["PATH"] = path + os.pathsep + os.environ["PATH"]
    elif os.path.isdir(path):
        os.environ["PATH"] = os.path.join(path, "bin") + os.pathsep + os.environ["PATH"]
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS
    return ContinuousStructure, openEMS


def model_bounds(project: Project, variables) -> tuple[np.ndarray, np.ndarray]:
    lo = np.full(3, np.inf)
    hi = np.full(3, -np.inf)
    for obj in project.objects:
        if obj.kind == "cad" and obj.material == "vacuum":
            continue  # transparent construction sheets from imported designs
        a, b = obj.bounds(variables, project)
        lo = np.minimum(lo, a)
        hi = np.maximum(hi, b)
    for port in project.ports:
        a, b = port.corners(variables)
        lo = np.minimum(lo, np.minimum(a, b))
        hi = np.maximum(hi, np.maximum(a, b))
    if not np.all(np.isfinite(lo)):
        lo, hi = -np.ones(3), np.ones(3)
    return lo, hi


def region_bounds(project: Project, variables) -> tuple[np.ndarray, np.ndarray]:
    """Simulation domain: model bounding box plus padding (lambda/4 at f_start when 'auto')."""
    if project.region.explicit:
        e = project.region.explicit
        return np.array([e[0], e[2], e[4]], float), np.array([e[1], e[3], e[5]], float)
    lo, hi = model_bounds(project, variables)
    pad_expr = project.region.padding.strip()
    if pad_expr.lower() == "auto":
        f_lo = max(float(project.setup.f_start), 0.25 * float(project.setup.f_stop)) * 1e9
        pad = C0 / f_lo / 4 / project.unit_scale()
    else:
        pad = evaluate(pad_expr, variables)
    return lo - pad, hi + pad


def _merge_close(lines, hard, tol, protected=()):
    """Collapse mesh lines closer than `tol`.

    A cluster keeps its 'protected' line (one face of a small gap between two conductors), otherwise its
    'hard' line (sheet/port plane), otherwise its midpoint. Two protected lines never share a cluster,
    so a gap narrower than `tol` keeps a mesh line on each side and stays open.
    """
    lines = np.unique(np.round(np.asarray(lines, dtype=float), 9))
    hard = np.unique(np.round(np.asarray(hard, dtype=float), 9))
    prot = np.unique(np.round(np.asarray(protected, dtype=float), 9))

    def is_prot(x):
        return len(prot) and np.any(np.abs(prot - x) < 1e-9)

    out, cluster = [], [lines[0]] if len(lines) else []
    for x in lines[1:]:
        if x - cluster[-1] < tol and not (is_prot(x) and any(is_prot(c) for c in cluster)):
            cluster.append(x)
            continue
        out += _representative(cluster, hard, prot)
        cluster = [x]
    if cluster:
        out += _representative(cluster, hard, prot)
    return np.array(out)


def _representative(cluster, hard, prot=()):
    """Lines to keep for one cluster: its protected lines (gap faces) and hard lines (sheet/port planes),
    all of them, since averaging two distinct planes would move both off their geometry (a thin
    substrate under a coarse grid used to lose its ground and top planes this way); else the midpoint."""
    keep = sorted({float(x) for x in cluster
                   for group in (prot, hard) if len(group) and np.any(np.abs(group - x) < 1e-9)})
    if keep:
        return keep
    return [float(cluster[0] if len(cluster) == 1 else np.mean([cluster[0], cluster[-1]]))]


def _collapse(values, feat):
    """Merge values closer than `feat` (chain clustering) into their cluster's midpoint."""
    v = np.unique(np.round(np.asarray(values, dtype=float), 9))
    if feat <= 0 or len(v) < 2:
        return list(v)
    out, cl = [], [v[0]]
    for x in v[1:]:
        if x - cl[-1] < feat:
            cl.append(x)
        else:
            out.append(0.5 * (cl[0] + cl[-1])); cl = [x]
    out.append(0.5 * (cl[0] + cl[-1]))
    return out


def _smooth(lines, max_res, hard=(), tol=0.0, protected=()):
    """Smooth mesh lines so no cell exceeds max_res and neighbours grow by <= 1.4x."""
    from CSXCAD.SmoothMeshLines import SmoothMeshLines
    lines = np.unique(np.round(np.asarray(list(lines) + list(protected), dtype=float), 9))
    if tol > 0:
        lines = _merge_close(lines, list(hard) + [lines[0], lines[-1]], tol, protected)
    return SmoothMeshLines(lines, max_res, 1.4)


# --------------------------------------------------------------------------- topology guards

def _conductor_boxes(project: Project, variables):
    """[(label, lo, hi)] for every conductor and every port terminal region."""
    out = []
    for obj in project.objects:
        if not project.material(obj.material).conductor:
            continue
        lo, hi = obj.bounds(variables, project)
        out.append((obj.name, np.array(lo, float), np.array(hi, float)))
    for port in project.ports:
        a, b = port.corners(variables)
        out.append((f"port {port.name}", np.minimum(a, b), np.maximum(a, b)))
    return out


def small_gaps(project: Project, variables, tol: float) -> list[dict]:
    """Pairs of conductors / port terminals that are closer than `tol` without touching.

    Merging mesh lines across such a gap would short it (or, for a port, connect a terminal that is
    floating in the source model). For each pair the axis with the largest separation is reported; one
    protected line on each side of the gap along that axis is enough to keep it open.
    """
    boxes = _conductor_boxes(project, variables)
    out = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            (na, la, ha), (nb, lb, hb) = boxes[i], boxes[j]
            if na.startswith("port ") and nb.startswith("port "):
                continue
            gap = np.maximum(np.maximum(lb - ha, la - hb), 0.0)
            d = float(np.linalg.norm(gap))
            if d >= tol:
                continue
            # bounding boxes are only a coarse filter; confirm with the exact distance between the shapes
            exact = _exact_distance(project, variables, na, nb)
            if exact is not None:
                d = exact
            if not 0 < d < tol:
                continue
            k = int(np.argmax(gap))
            if lb[k] - ha[k] > 0:
                faces = (float(ha[k]), float(lb[k]))
            else:
                faces = (float(hb[k]), float(la[k]))
            out.append({"a": na, "b": nb, "axis": k, "faces": faces, "distance": d})
    return out


def _occ_shape(project: Project, variables, label):
    """OpenCascade shape for an object or a port region (native primitives use their box)."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt
    if label.startswith("port "):
        port = next(q for q in project.ports if q.name == label[5:])
        a, b = (np.array(v, float) for v in port.corners(variables))
        lo, hi = np.minimum(a, b), np.maximum(a, b)
    else:
        obj = next(o for o in project.objects if o.name == label)
        if obj.kind == "cad":
            from . import cad
            return cad.shape(project, obj)
        lo, hi = (np.array(v, float) for v in obj.bounds(variables, project))
    ext = hi - lo
    if np.count_nonzero(ext > 1e-12) == 3:
        return BRepPrimAPI_MakeBox(gp_Pnt(*lo), gp_Pnt(*hi)).Shape()
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeFace
    flat = [i for i in range(3) if ext[i] <= 1e-12]
    if len(flat) == 1:
        u, w = [i for i in range(3) if i != flat[0]]
        c = [lo.copy(), lo.copy(), lo.copy(), lo.copy()]
        c[1][u] = hi[u]
        c[2][u], c[2][w] = hi[u], hi[w]
        c[3][w] = hi[w]
        poly = BRepBuilderAPI_MakePolygon(*[gp_Pnt(*q) for q in c], True)
        return BRepBuilderAPI_MakeFace(poly.Wire(), True).Face()
    return BRepBuilderAPI_MakeEdge(gp_Pnt(*lo), gp_Pnt(*hi)).Edge()


def _exact_distance(project: Project, variables, la, lb):
    try:
        from .aedt.geometry import dist
        return float(dist(_occ_shape(project, variables, la), _occ_shape(project, variables, lb)))
    except Exception:  # noqa: BLE001 - fall back to the bounding-box estimate
        return None


def port_connectivity(project: Project, variables) -> list[dict]:
    """Distance from each port terminal to the nearest conductor (0 = connected).

    Uses the exact B-rep distance for imported objects and the bounding box for native primitives.
    """
    conductors = [o for o in project.objects if project.material(o.material).conductor]
    lo_m, hi_m = model_bounds(project, variables)
    touch = 1e-6 * float(np.linalg.norm(hi_m - lo_m)) + 1e-9
    out = []
    for port in project.ports:
        a, b = (np.array(v, float) for v in port.corners(variables))
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        ax = AXES.index(port.axis)
        centre = 0.5 * (lo + hi)
        for label, coord in (("start", lo[ax]), ("end", hi[ax])):
            pt = centre.copy()
            pt[ax] = coord
            best, best_d = None, float("inf")
            for obj in conductors:
                d = _point_distance(project, obj, pt, variables)
                if d < best_d:
                    best, best_d = obj.name, d
            out.append({"port": port.name, "terminal": label, "point": pt.tolist(), "nearest": best,
                        "distance": best_d, "connected": best_d <= touch})
    return out


def _point_distance(project, obj, pt, variables) -> float:
    if obj.kind == "cad":
        from . import cad
        from .aedt.geometry import dist, vertex
        return dist(vertex(pt), cad.shape(project, obj))
    lo, hi = (np.array(v, float) for v in obj.bounds(variables, project))
    return float(np.linalg.norm(np.maximum(np.maximum(lo - pt, pt - hi), 0.0)))


def connectivity_warnings(project: Project, variables) -> list[str]:
    msgs = []
    for c in port_connectivity(project, variables):
        if not c["connected"]:
            msgs.append(f"Port '{c['port']}' {c['terminal']} terminal does not touch a conductor: nearest is "
                        f"'{c['nearest']}' at {c['distance'] * 1000:.1f} um ({project.units} model). The port "
                        f"couples only capacitively; check this matches the source design")
    return msgs


def build_mesh(project: Project, variables, rlo, rhi, scale: float = 1.0) -> list[np.ndarray]:
    """Mesh lines per axis. `scale` > 1 refines everything proportionally (adaptive passes)."""
    s = project.setup
    f_max = float(s.f_stop) * 1e9
    eps_max = max([project.material(o.material).epsr * project.material(o.material).mur
                   for o in project.objects if not project.material(o.material).conductor] + [1.0])
    eps_max = min(eps_max, 100.0)
    lam_min = C0 / f_max / np.sqrt(eps_max) / project.unit_scale()
    max_res = lam_min / (s.cells_per_wavelength * scale)
    n_thin = int(np.ceil(4 * scale))   # cells across thin dielectric layers

    lines = [[rlo[i], rhi[i]] for i in range(3)]
    hard = [[], [], []]   # planes that must stay exactly on a mesh line (sheets, ports)
    min_edge = max_res    # finest edge cell requested (keeps line merging from undoing it)
    min_layer = max_res   # finest thin-layer cell requested (likewise)
    for obj in project.objects:
        mat = project.material(obj.material)
        if obj.kind == "cad" and obj.material in ("vacuum", "air") and not mat.conductor:
            continue
        a, b = obj.bounds(variables, project)
        for i in range(3):
            lines[i] += [a[i], b[i]]
            thick = b[i] - a[i]
            if thick < 1e-9 * max(1.0, abs(a[i])):
                hard[i].append(a[i])
            # resolve thin dielectric layers (substrates) with at least four cells (more when refined)
            if 0 < thick < n_thin * max_res and not mat.conductor:
                lines[i] += list(np.linspace(a[i], b[i], n_thin + 1))
                min_layer = min(min_layer, thick / n_thin)
            if obj.kind in ("cylinder", "sphere") and thick > 0:
                lines[i] += list(np.linspace(a[i], b[i], 7))
        # Conductor sheets: FDTD error concentrates at their edges (field singularity). Validation on a
        # 50 ohm microstrip showed eps_eff +8.6 % with edge cells of res/4, +0.5 % with res/8, and the
        # wavelength-based cell can be wider than the strip itself, so size the edge cells by the
        # sheet's own width as well as by the wavelength.
        if s.edge_refinement and mat.conductor:
            ext = np.asarray(b, float) - np.asarray(a, float)
            flat = ext < 1e-9 * max(1.0, float(np.abs(np.asarray(a, float)).max()))
            if flat.sum() == 1:
                width = float(ext[~flat].min())
                edge = min(max_res / 4, width / (10 * scale))
                min_edge = min(min_edge, edge)
                for i in np.flatnonzero(~flat):
                    lines[i] += [a[i] - edge, a[i] + edge, b[i] - edge, b[i] + edge]
        if obj.kind == "cad" and mat.conductor:
            from . import cad
            for i, coords in enumerate(cad.vertex_coords(project, obj)):
                lines[i] += list(coords)
    for port in project.ports:
        a, b = port.corners(variables)
        for i in range(3):
            lines[i] += [a[i], b[i]]
            hard[i] += [a[i], b[i]]
            # A port sheet with width (imported HFSS lumped/wave ports) needs a line through its centre:
            # openEMS puts the voltage probe there and fails ("only a 1D/line integration is allowed")
            # if the probe falls between two mesh lines.
            if AXES[i] != port.axis.upper() and abs(b[i] - a[i]) > 0:
                lines[i].append(0.5 * (a[i] + b[i]))
                hard[i].append(0.5 * (a[i] + b[i]))
    tol = min(max_res / 10, min_edge / 2, min_layer / 2)
    # never merge across a gap between two conductors (or a port terminal and a conductor):
    # that would change the circuit topology, e.g. connect a port that is floating in the source design
    protected = [[], [], []]
    for g in small_gaps(project, variables, max(tol, max_res / 10)):
        protected[g["axis"]] += list(g["faces"])
    feat = float(getattr(s, "min_feature", 0.0) or 0.0)
    if feat > 0:
        # Unresolved features (speed option): fixed planes and protected gap faces closer than `feat`
        # become one line at their midpoint; sheets are thickened to >= feat so they still cover it.
        tol = max(tol, feat)
        for i in range(3):
            hard[i] = _collapse(hard[i], feat)
            protected[i] = _collapse(protected[i], feat)
            fixed = np.asarray(hard[i] + protected[i])
            if len(fixed):
                lines[i] = [x for x in lines[i] if np.min(np.abs(fixed - x)) >= feat] + list(fixed)
    lines = [np.clip(lines[i], rlo[i], rhi[i]) for i in range(3)]
    return [_smooth(lines[i], max_res, hard[i], tol, protected[i]) for i in range(3)], max_res


def _slab_for(project, mesh, max_res):
    """Thickness given to zero-thickness sheets (polyhedra need volume). max_res/5 covers planes that
    line merging moved by up to max_res/10; with min_feature set, at least that feature size."""
    feat = float(getattr(project.setup, "min_feature", 0.0) or 0.0)
    return max(max_res / 5, 1.01 * feat)


def build_simulation(project: Project, sim_path: str | None = None, scale: float = 1.0):
    """Return (FDTD, CSX, ports, nf2ff, info dict, mesh). `scale` refines the mesh (adaptive passes)."""
    ContinuousStructure, openEMS = _import_openems()
    variables = project.resolved_variables()
    s = project.setup
    f_start, f_stop = float(s.f_start) * 1e9, float(s.f_stop) * 1e9
    f0 = 0.5 * (f_start + f_stop)
    fc = 0.5 * (f_stop - f_start)

    fdtd = openEMS(NrTS=int(s.max_timesteps), EndCriteria=10 ** (s.end_criteria_db / 10))
    fdtd.SetGaussExcite(f0, fc)
    fdtd.SetBoundaryCond(list(project.region.boundaries))
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(project.unit_scale())

    rlo, rhi = region_bounds(project, variables)
    mesh, max_res = build_mesh(project, variables, rlo, rhi, scale)
    for i, ax in enumerate("xyz"):
        grid.SetLines(ax, mesh[i])

    # Timestep budget from the mesh. The FDTD timestep is set by the smallest cells (Courant limit), so
    # fine features (e.g. 23 um copper offsets in imported designs) can make a fixed budget too short to
    # even finish the excitation pulse. Allow 20 excitation lengths; openEMS stops earlier once the energy
    # has decayed to the end criterion.
    unit = project.unit_scale()
    dmin = [float(np.diff(m).min()) * unit for m in mesh]
    dt_est = 1.0 / (C0 * np.sqrt(sum(1.0 / d ** 2 for d in dmin)))
    t_exc = 2 * 9 / (2 * np.pi * fc)                # openEMS Gauss pulse length
    nrts = max(int(s.max_timesteps), int(20 * t_exc / dt_est))
    fdtd.SetNumberOfTimeSteps(nrts)

    # Materials: one CSX property per (material) and one primitive per object.
    props = {}
    for idx, obj in enumerate(project.objects):
        mat = project.material(obj.material)
        if obj.material in ("vacuum", "air") and not mat.conductor:
            continue  # background medium, nothing to add
        if obj.material not in props:
            if mat.conductor:
                props[obj.material] = csx.AddMetal(_safe(obj.material))
            else:
                kappa = 2 * np.pi * f0 * EPS0 * mat.epsr * mat.tand
                props[obj.material] = csx.AddMaterial(_safe(obj.material), epsilon=mat.epsr,
                                                      mue=mat.mur, kappa=kappa)
        prop = props[obj.material]
        prio = idx + (1000 if mat.conductor else 0)  # metal always wins over dielectric
        v = obj.values(variables)
        c = [v.get("x", 0.0), v.get("y", 0.0), v.get("z", 0.0)]
        a = AXES.index(obj.axis)
        if obj.kind in ("box", "rectangle"):
            lo, hi = obj.bounds(variables)
            prop.AddBox(lo, hi, priority=prio)
        elif obj.kind == "cylinder":
            stop = list(c)
            stop[a] += v["height"]
            prop.AddCylinder(c, stop, v["radius"], priority=prio)
        elif obj.kind == "sphere":
            prop.AddSphere(c, v["radius"], priority=prio)
        elif obj.kind == "cad":
            from . import cad
            # sheets get a thickness of two snapping tolerances so the mesh plane they sit on is inside
            pts, tris = cad.mesh(project, obj, slab=_slab_for(project, mesh, max_res))
            if len(tris) == 0:
                continue
            poly = prop.AddPolyhedron(priority=prio)
            for q in pts:
                poly.AddVertex(float(q[0]), float(q[1]), float(q[2]))
            for t in tris:
                poly.AddFace([int(t[0]), int(t[1]), int(t[2])])

    ports = []
    for n, port in enumerate(project.ports, start=1):
        a, b = port.corners(variables)
        r = evaluate(port.impedance, variables)
        ports.append(fdtd.AddLumpedPort(n, r, a, b, port.axis.lower(), excite=1 if port.excited else 0))

    nf2ff = None
    if s.farfield:
        nf2ff = fdtd.CreateNF2FFBox()

    info = {
        "cells": int(np.prod([len(m) for m in mesh])),
        "mesh_lines": [len(m) for m in mesh],
        "max_res": max_res,
        "min_cell": min(float(np.diff(m).min()) for m in mesh if len(m) > 1),
        "region": (rlo.tolist(), rhi.tolist()),
        "f0": f0, "fc": fc, "scale": scale, "nrts": nrts, "dt_est": dt_est,
    }
    return fdtd, csx, ports, nf2ff, info, mesh


def _safe(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in name)


def solve(project_file: str, out_dir: str) -> str:
    project = Project.load(project_file)
    out_dir = os.path.abspath(out_dir)
    sim_path = os.path.join(out_dir, "sim")
    os.makedirs(sim_path, exist_ok=True)
    s = project.setup

    if not project.ports:
        raise RuntimeError("No excitation defined: add at least one lumped port.")
    if not any(p.excited for p in project.ports):
        raise RuntimeError("No port is excited.")

    variables = project.resolved_variables()
    for msg in connectivity_warnings(project, variables):
        print(f"[FieldForge] WARNING: {msg}", flush=True)

    def run_pass(scale: float):
        fdtd, csx, ports, nf2ff, info, mesh = build_simulation(project, sim_path, scale)
        for g in small_gaps(project, variables, info["max_res"] / 10):
            print(f"[FieldForge] Kept a {g['distance'] * 1000:.1f} um gap between '{g['a']}' and '{g['b']}' "
                  f"open (separate mesh lines along {'xyz'[g['axis']]})", flush=True)
        print(f"[FieldForge] Mesh: {info['mesh_lines'][0]} x {info['mesh_lines'][1]} x "
              f"{info['mesh_lines'][2]} = {info['cells']:,} cells "
              f"(cells {info['min_cell']:.3g}-{info['max_res']:.3g} {project.units})", flush=True)
        note = "" if info["nrts"] == int(s.max_timesteps) else f" (raised from {s.max_timesteps} for the mesh's small cells)"
        print(f"[FieldForge] Max timesteps {info['nrts']:,}{note}, timestep ~{info['dt_est'] * 1e12:.3g} ps, "
              f"end criteria {s.end_criteria_db} dB", flush=True)
        ret = fdtd.Run(sim_path, cleanup=True, verbose=1)
        if ret:
            raise RuntimeError(f"openEMS returned code {ret}: usually the timestep limit ({info['nrts']:,}) was reached "
                               f"before the field energy decayed to {s.end_criteria_db} dB, so the S-parameters would "
                               f"be truncated and wrong. Raise 'Max time steps' or relax the end criterion.")
        print("[FieldForge] Post-processing ports...", flush=True)
        for p in ports:
            p.CalcPort(sim_path, f)
        inc = ports[exc].uf_inc
        # Column `exc` of the scattering matrix: outgoing wave at every port / incident wave at the excited one.
        S = np.array([p.uf_ref / inc for p in ports])
        return dict(S=S, ports=ports, nf2ff=nf2ff, mesh=mesh, n_cells=info["cells"])

    f = np.linspace(float(s.f_start) * 1e9, float(s.f_stop) * 1e9, int(s.n_points))
    exc = next(i for i, p in enumerate(project.ports) if p.excited)
    history = None
    if s.adaptive:
        from .adaptive import adaptive_mesh

        def estimate(scale):
            variables = project.resolved_variables()
            rlo, rhi = region_bounds(project, variables)
            return int(np.prod([len(m) for m in build_mesh(project, variables, rlo, rhi, scale)[0]]))

        print(f"[FieldForge] Adaptive passes: max {s.max_passes}, Maximum Delta S {s.max_delta_s:g}, "
              f"min converged passes {s.min_converged_passes}", flush=True)
        ad = adaptive_mesh(lambda sc: run_pass(sc), estimate, start_scale=1.0,
                           refine_factor=float(s.refine_factor), max_delta=float(s.max_delta_s),
                           min_converged_passes=int(s.min_converged_passes), max_passes=int(s.max_passes),
                           max_cells=int(s.max_cells),
                           log=lambda m: print(m.replace("[adaptive]", "[FieldForge] Adaptive"), flush=True))
        if not ad.passes:
            raise RuntimeError(f"Adaptive meshing could not start: {ad.reason}")
        print(f"[FieldForge] Adaptive {'converged' if ad.converged else 'NOT converged'}: {ad.reason}", flush=True)
        out = ad.final.output
        history = np.array([(p.index, p.scale, p.cells, np.nan if p.delta_s is None else p.delta_s, p.runtime_s)
                            for p in ad.passes])
    else:
        out = run_pass(1.0)
    S, ports, nf2ff, mesh = out["S"], out["ports"], out["nf2ff"], out["mesh"]
    if np.all(np.abs(1 - np.abs(S[exc])) < 1e-6):
        # no current flowed: the result is the trivial open circuit, not a measurement
        print(f"[FieldForge] WARNING: port '{project.ports[exc].name}' carried no current (|S11| = 1 at every "
              f"frequency). The port is not driving the structure; check that its terminals and mesh lines "
              f"are placed correctly. These results are not meaningful.", flush=True)
    Zin = np.array([p.uf_tot / p.if_tot for p in ports])
    P_acc = np.array([np.real(p.P_acc) for p in ports])

    results = dict(freq=f, S=S, Zin=Zin, P_acc=P_acc, excited=exc,
                   port_names=np.array([p.name for p in project.ports]),
                   mesh_x=mesh[0], mesh_y=mesh[1], mesh_z=mesh[2],
                   unit=project.unit_scale())
    if history is not None:
        # columns: pass, mesh scale, cells, max|dS| vs previous pass (NaN for the first), runtime (s)
        results.update(adaptive_history=history, adaptive_converged=bool(ad.converged))

    if nf2ff is not None:
        ff_f = float(s.farfield_freq) * 1e9
        step = float(s.farfield_step_deg)
        theta = np.arange(0.0, 180.0 + step / 2, step)
        phi = np.arange(0.0, 360.0 + step / 2, step)
        print(f"[FieldForge] Computing far field at {ff_f / 1e9:g} GHz...", flush=True)
        # the phase centre must lie inside the near-field recording box: use the model centre
        mlo, mhi = model_bounds(project, project.resolved_variables())
        centre = (0.5 * (mlo + mhi) * project.unit_scale()).tolist()
        res = nf2ff.CalcNF2FF(sim_path, ff_f, theta, phi, center=centre, read_cached=False)
        e_norm = res.E_norm[0]
        dmax = float(np.real(res.Dmax[0]))
        D = (e_norm / e_norm.max()) ** 2 * dmax
        P_in = float(np.interp(ff_f, f, P_acc[exc]))
        results.update(ff_freq=ff_f, ff_theta=theta, ff_phi=phi, ff_D=D, ff_Dmax=dmax,
                       ff_Prad=float(np.real(res.Prad[0])), ff_Pin=P_in,
                       ff_Etheta=np.abs(res.E_theta[0]), ff_Ephi=np.abs(res.E_phi[0]))
        s11_ff = abs(np.interp(ff_f, f, np.abs(S[exc])))
        accepted = 1 - s11_ff ** 2
        results.update(ff_accepted=accepted)
        if accepted < MIN_ACCEPTED:
            eff_txt = (f"radiation efficiency not meaningful (port accepts only {100 * accepted:.2f} % of the "
                       f"incident power at {ff_f / 1e9:g} GHz)")
        else:
            eff_txt = f"radiation efficiency {100 * float(np.real(res.Prad[0])) / max(P_in, 1e-30):.1f} %"
        print(f"[FieldForge] Max directivity {10 * np.log10(dmax):.2f} dBi, {eff_txt}", flush=True)

    out_file = os.path.join(out_dir, "results.npz")
    np.savez(out_file, **results)
    print(f"[FieldForge] Results written to {out_file}", flush=True)
    return out_file


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    try:
        solve(sys.argv[1], sys.argv[2])
    except Exception as exc:  # report cleanly to the GUI log
        import traceback
        traceback.print_exc()
        print(f"[FieldForge] ERROR: {exc}", flush=True)
        sys.exit(1)
