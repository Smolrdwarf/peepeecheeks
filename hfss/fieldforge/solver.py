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
    return os.environ.get("OPENEMS_INSTALL_PATH",
                          os.path.join(os.path.expanduser("~"), "tools", "openEMS"))


def _import_openems():
    path = openems_install_path()
    os.environ["OPENEMS_INSTALL_PATH"] = path
    if os.name == "nt" and os.path.isdir(path):
        os.add_dll_directory(path)
        os.environ["PATH"] = path + os.pathsep + os.environ["PATH"]
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
        out.append(_representative(cluster, hard, tol))
        cluster = [x]
    if cluster:
        out.append(_representative(cluster, hard, tol))
    return np.array(out)


def _representative(cluster, hard, tol):
    h = [x for x in cluster if np.any(np.abs(hard - x) < 1e-9)]
    if h:
        return float(np.mean(h))
    return float(cluster[0] if len(cluster) == 1 else np.mean([cluster[0], cluster[-1]]))


def _smooth(lines, max_res, hard=(), tol=0.0):
    """Smooth mesh lines so no cell exceeds max_res and neighbours grow by <= 1.4x."""
    from CSXCAD.SmoothMeshLines import SmoothMeshLines
    lines = np.unique(np.round(np.asarray(lines, dtype=float), 9))
    if tol > 0:
        lines = _merge_close(lines, list(hard) + [lines[0], lines[-1]], tol)
    return SmoothMeshLines(lines, max_res, 1.4)


def build_mesh(project: Project, variables, rlo, rhi) -> list[np.ndarray]:
    s = project.setup
    f_max = float(s.f_stop) * 1e9
    eps_max = max([project.material(o.material).epsr * project.material(o.material).mur
                   for o in project.objects if not project.material(o.material).conductor] + [1.0])
    eps_max = min(eps_max, 100.0)
    lam_min = C0 / f_max / np.sqrt(eps_max) / project.unit_scale()
    max_res = lam_min / s.cells_per_wavelength

    lines = [[rlo[i], rhi[i]] for i in range(3)]
    hard = [[], [], []]   # planes that must stay exactly on a mesh line (sheets, ports)
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
            # resolve thin dielectric layers (substrates) with at least four cells
            if 0 < thick < 4 * max_res and not mat.conductor:
                lines[i] += list(np.linspace(a[i], b[i], 5))
            if obj.kind in ("cylinder", "sphere") and thick > 0:
                lines[i] += list(np.linspace(a[i], b[i], 7))
        if obj.kind == "cad" and mat.conductor:
            from . import cad
            for i, coords in enumerate(cad.vertex_coords(project, obj)):
                lines[i] += list(coords)
    for port in project.ports:
        a, b = port.corners(variables)
        for i in range(3):
            lines[i] += [a[i], b[i]]
            hard[i] += [a[i], b[i]]
    tol = max_res / 10
    return [_smooth(lines[i], max_res, hard[i], tol) for i in range(3)], max_res


def build_simulation(project: Project, sim_path: str | None = None):
    """Return (FDTD, CSX, ports, nf2ff, info dict)."""
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
    mesh, max_res = build_mesh(project, variables, rlo, rhi)
    for i, ax in enumerate("xyz"):
        grid.SetLines(ax, mesh[i])

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
            pts, tris = cad.mesh(project, obj, slab=max_res / 5)
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
        "region": (rlo.tolist(), rhi.tolist()),
        "f0": f0, "fc": fc,
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

    fdtd, csx, ports, nf2ff, info, mesh = build_simulation(project, sim_path)
    print(f"[FieldForge] Mesh: {info['mesh_lines'][0]} x {info['mesh_lines'][1]} x "
          f"{info['mesh_lines'][2]} = {info['cells']:,} cells "
          f"(max cell {info['max_res']:.3g} {project.units})", flush=True)
    print(f"[FieldForge] NrTS={s.max_timesteps}  end criteria {s.end_criteria_db} dB", flush=True)

    ret = fdtd.Run(sim_path, cleanup=True, verbose=1)
    if ret:
        raise RuntimeError(f"openEMS setup failed (code {ret})")

    print("[FieldForge] Post-processing ports...", flush=True)
    f = np.linspace(float(s.f_start) * 1e9, float(s.f_stop) * 1e9, int(s.n_points))
    for p in ports:
        p.CalcPort(sim_path, f)

    exc = next(i for i, p in enumerate(project.ports) if p.excited)
    inc = ports[exc].uf_inc
    # Column `exc` of the scattering matrix: outgoing wave at every port / incident wave at the excited one.
    S = np.array([p.uf_ref / inc for p in ports])
    Zin = np.array([p.uf_tot / p.if_tot for p in ports])
    P_acc = np.array([np.real(p.P_acc) for p in ports])

    results = dict(freq=f, S=S, Zin=Zin, P_acc=P_acc, excited=exc,
                   port_names=np.array([p.name for p in project.ports]),
                   mesh_x=mesh[0], mesh_y=mesh[1], mesh_z=mesh[2],
                   unit=project.unit_scale())

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
