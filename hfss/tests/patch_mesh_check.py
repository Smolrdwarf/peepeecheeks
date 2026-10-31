"""Patch-antenna example: S11 resonance with edge refinement off / on / adaptive.  ~10 min, needs openEMS."""
import os
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fieldforge import examples              # noqa: E402
from fieldforge.solver import solve          # noqa: E402

work = os.path.join(tempfile.gettempdir(), "fieldforge_patch_check")
os.makedirs(work, exist_ok=True)
for name, kw in [("edge_off", dict(edge_refinement=False)), ("edge_on", dict(edge_refinement=True)),
                 ("adaptive", dict(edge_refinement=True, adaptive=True, max_cells=4_000_000))]:
    p = examples.patch_antenna()
    p.setup.farfield = False
    for k, v in kw.items():
        setattr(p.setup, k, v)
    pf = os.path.join(work, f"{name}.ffproj")
    p.save(pf)
    t0 = time.time()
    d = np.load(solve(pf, os.path.join(work, name)))
    s11 = 20 * np.log10(np.abs(d["S"][d["excited"]]))
    i = int(np.argmin(s11))
    cells = len(d["mesh_x"]) * len(d["mesh_y"]) * len(d["mesh_z"])
    extra = ""
    if "adaptive_history" in d:
        extra = "  passes: " + "; ".join(f"x{r[1]:.2f} {int(r[2]):,} dS={r[3]:.3f}" for r in d["adaptive_history"])
    print(f">>> {name}: f_res {d['freq'][i] / 1e9:.4f} GHz, S11 min {s11[i]:.1f} dB, {cells:,} cells, "
          f"{time.time() - t0:.0f} s{extra}", flush=True)
