"""Run one openEMS variant of the imported ring antenna with a single mechanism knob changed.

usage: python run_variant.py <project.ffproj> <out_dir> <name> [knob=value ...]

knobs (defaults reproduce the FieldForge baseline):
  cpw=20          cells per wavelength
  gap=none        port-to-ground gap in mm: keep separate mesh lines at the port plane and at
                  port+gap, and carve the ground plane back to x = port+gap (HFSS geometry: 0.0427)
  tand=1          scale factor on the dielectric loss tangent (0 = lossless)
  pad=none        region padding in mm on every side (overrides the imported HFSS region)
  ec=-40          end criterion (energy decay, dB)
  slab=0.2        sheet slab thickness as a fraction of the max cell (baseline max_res/5)
"""
import json
import os
import sys
import time

import numpy as np

HFSS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, HFSS_DIR)

from fieldforge import solver  # noqa: E402
from fieldforge.model import Project  # noqa: E402


def main():
    proj_file, out_dir, name = sys.argv[1:4]
    out_dir = os.path.abspath(out_dir)
    knobs = dict(a.split("=", 1) for a in sys.argv[4:])
    p = Project.load(proj_file)
    p.setup.farfield = False
    p.setup.cells_per_wavelength = int(knobs.get("cpw", 20))
    p.setup.end_criteria_db = float(knobs.get("ec", -40))
    if "tand" in knobs:
        for m in p.materials.values():
            m.tand *= float(knobs["tand"])
    if knobs.get("pad", "none") != "none":
        p.region.explicit = None
        p.region.padding = knobs["pad"]

    slab_frac = float(knobs.get("slab", 0.2))
    from fieldforge import cad
    real_mesh = cad.mesh

    def mesh_with_slab(project, obj, slab=0.0):
        return real_mesh(project, obj, slab=slab * slab_frac / 0.2 if slab else 0.0)
    cad.mesh = mesh_with_slab

    sim = os.path.abspath(os.path.join(out_dir, name, "sim"))
    os.makedirs(sim, exist_ok=True)
    fdtd, csx, ports, nf2ff, info, mesh = solver.build_simulation(p, sim)

    gap = knobs.get("gap", "none")
    if gap != "none":
        from CSXCAD.SmoothMeshLines import SmoothMeshLines
        g = float(gap)
        port = p.ports[0]
        xp = float(port.params["x1"])
        xg = xp + g
        x = np.asarray(mesh[0])
        x = np.concatenate([x[(x < xp - 1e-9) | (x > xg + 1e-9)], [xp, xg]])
        x = SmoothMeshLines(np.unique(np.round(x, 9)), info["max_res"], 1.4)
        csx.GetGrid().SetLines("x", x)
        mesh[0] = x
        # carve the ground plane back to x = xg with a high-priority vacuum box (ground sits at z = -0.03)
        carve = csx.AddMaterial("carve_vacuum", epsilon=1.0)
        e = p.region.explicit or [x[0], x[-1], mesh[1][0], mesh[1][-1], mesh[2][0], mesh[2][-1]]
        carve.AddBox([x[0] - 1, e[2] - 1, -0.03 - 0.25], [xg, e[3] + 1, -0.03 + 0.25], priority=100000)

    cells = int(np.prod([len(m) for m in mesh]))
    dmin = min(np.diff(m).min() for m in mesh)
    t0 = time.time()
    ret = fdtd.Run(sim, cleanup=True, verbose=0)
    if ret:
        raise RuntimeError(f"openEMS failed ({ret})")
    s = p.setup
    f = np.linspace(float(s.f_start) * 1e9, float(s.f_stop) * 1e9, int(s.n_points))
    ports[0].CalcPort(sim, f)
    S11 = ports[0].uf_ref / ports[0].uf_inc
    Z = ports[0].uf_tot / ports[0].if_tot
    np.savez(os.path.join(out_dir, name, "s11.npz"), f=f, S11=S11, Z=Z)
    meta = {"name": name, "knobs": knobs, "cells": cells, "min_cell_mm": float(dmin),
            "mesh_lines": [len(m) for m in mesh], "runtime_s": round(time.time() - t0, 1)}
    json.dump(meta, open(os.path.join(out_dir, name, "meta.json"), "w"), indent=1)
    k = int(np.argmin(np.abs(S11)))
    print(f"{name}: {cells:,} cells, min cell {dmin:.4f} mm, {meta['runtime_s']} s; "
          f"min |S11| {20 * np.log10(abs(S11[k])):.2f} dB at {f[k] / 1e9:.3f} GHz")


if __name__ == "__main__":
    main()
