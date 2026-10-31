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


def _merge_close(lines, hard, tol):
    """Collapse mesh lines closer than `tol`; a cluster keeps its 'hard' line (sheet/port plane) if any."""
    lines = np.unique(np.round(np.asarray(lines, dtype=float), 9))
    hard = np.unique(np.round(np.asarray(hard, dtype=float), 9))
    out, cluster = [], [lines[0]] if len(lines) else []
    for x in lines[1:]:
        if x - cluster[-1] < tol:
            cluster.append(x)
            continue
        out += _representative(cluster, hard, tol)
        cluster = [x]
    if cluster:
        out += _representative(cluster, hard, tol)
    return np.array(out)


def _representative(cluster, hard, tol):
    """Lines to keep for one cluster. Hard lines are distinct planes (sheets, port ends, region faces)
    that must all survive; averaging two of them would move both off their geometry."""
    h = [float(x) for x in cluster if np.any(np.abs(hard - x) < 1e-9)]
    if h:
        return h
    return [float(cluster[0] if len(cluster) == 1 else np.mean([cluster[0], cluster[-1]]))]


def _smooth(lines, max_res, hard=(), tol=0.0):
    """Smooth mesh lines so no cell exceeds max_res and neighbours grow by <= 1.4x."""
    from CSXCAD.SmoothMeshLines import SmoothMeshLines
    lines = np.unique(np.round(np.asarray(lines, dtype=float), 9))
    if tol > 0:
        lines = _merge_close(lines, list(hard) + [lines[0], lines[-1]], tol)
    return SmoothMeshLines(lines, max_res, 1.4)


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
            elif mat.conductor:
                # A conductor's extent must not be moved by line merging: snapping its edge onto a
                # neighbouring line leaves no metal on that line, so the conductor shrinks by up to a
                # cell and small gaps silently open or close. Example: an imported HFSS lumped port
                # 43 um short of the ground edge became a 1.4 mm gap and the port read as open.
                hard[i] += [a[i], b[i]]
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
    lines = [np.clip(lines[i], rlo[i], rhi[i]) for i in range(3)]
    return [_smooth(lines[i], max_res, hard[i], tol) for i in range(3)], max_res


def port_gaps(project: Project, variables=None) -> list[str]:
    """Warnings for lumped-port ends that do not touch any conductor.

    Uses conductor bounding boxes, which can only under-estimate a gap: a reported gap is real, while a
    port end inside a box (e.g. over a hole) is not flagged. A floating port end reads as an open circuit
    (or a small gap capacitance), which is usually a modelling mistake in the source design.
    """
    variables = project.resolved_variables() if variables is None else variables
    boxes = []
    for obj in project.objects:
        if project.material(obj.material).conductor:
            a, b = obj.bounds(variables, project)
            boxes.append((obj.name, np.asarray(a, float), np.asarray(b, float)))
    out = []
    for port in project.ports:
        a, b = (np.asarray(c, float) for c in port.corners(variables))
        ax = AXES.index(port.axis.upper())
        for end in (a[ax], b[ax]):
            pt = 0.5 * (a + b)
            pt[ax] = end
            if not boxes:
                continue
            dists = [(float(np.linalg.norm(np.maximum(0, np.maximum(lo - pt, pt - hi)))), name) for name, lo, hi in boxes]
            d, name = min(dists)
            if d > 1e-9 * max(1.0, float(np.abs(pt).max())):
                out.append(f"Port '{port.name}': the end at {port.axis.lower()} = {end:g} {project.units} does not touch "
                           f"a conductor (nearest: {name}, {d * project.unit_scale() * 1e6:.0f} um away). The port is "
                           f"floating there and will behave like an open circuit through a gap capacitance.")
    return out


def sheet_slab(mesh, max_res) -> float:
    """Thickness given to zero-thickness sheets so their mesh plane lies inside them. Kept below the
    smallest cell so a sheet never reaches the neighbouring line (which would close real gaps)."""
    min_cell = min(float(np.diff(m).min()) for m in mesh if len(m) > 1)
    return min(max_res / 5, 0.9 * min_cell)


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
            pts, tris = cad.mesh(project, obj, slab=sheet_slab(mesh, max_res))
            if len(tris) == 0:
                continue
            if mat.conductor:
                # CSXCAD treats points exactly on a polyhedron's surface as outside, so a conductor whose
                # edge sits on a mesh line (as every conductor extent does, see build_mesh) would have no
                # metal on that line. Grow it outward along each axis by 0.1 % of the smallest cell so
                # its bounding planes are inside. (Boxes/rectangles are inclusive already.)
                pts = np.asarray(pts, float)
                lo_p, hi_p = pts.min(axis=0), pts.max(axis=0)
                c, half = 0.5 * (lo_p + hi_p), 0.5 * (hi_p - lo_p)
                grow = 1e-3 * min(float(np.diff(m).min()) for m in mesh if len(m) > 1)
                k = np.where(half > grow, 1.0 + grow / np.maximum(half, 1e-30), 1.0)
                pts = c + (pts - c) * k
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

    for w in port_gaps(project):
        print(f"[FieldForge] WARNING: {w}", flush=True)

    def run_pass(scale: float):
        fdtd, csx, ports, nf2ff, info, mesh = build_simulation(project, sim_path, scale)
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
        res = nf2ff.CalcNF2FF(sim_path, ff_f, theta, phi, center=[0, 0, 0], read_cached=False)
        e_norm = res.E_norm[0]
        dmax = float(np.real(res.Dmax[0]))
        D = (e_norm / e_norm.max()) ** 2 * dmax
        P_in = float(np.interp(ff_f, f, P_acc[exc]))
        results.update(ff_freq=ff_f, ff_theta=theta, ff_phi=phi, ff_D=D, ff_Dmax=dmax,
                       ff_Prad=float(np.real(res.Prad[0])), ff_Pin=P_in,
                       ff_Etheta=np.abs(res.E_theta[0]), ff_Ephi=np.abs(res.E_phi[0]))
        print(f"[FieldForge] Max directivity {10 * np.log10(dmax):.2f} dBi, "
              f"radiation efficiency {100 * float(np.real(res.Prad[0])) / max(P_in, 1e-30):.1f} %",
              flush=True)

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
