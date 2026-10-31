"""Does port current flow?  python port_diag.py <project> <outdir> <variant>"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
from fieldforge import solver  # noqa: E402
from fieldforge.model import Project  # noqa: E402

proj, out, variant = sys.argv[1:4]
out = os.path.abspath(out)
p = Project.load(proj)
p.setup.farfield = False
p.setup.max_timesteps = 3000
p.setup.cells_per_wavelength = 10
if variant == "E3":
    solver.small_gaps = lambda *a, **k: []          # old behaviour: no protected lines
sim = os.path.join(out, variant)
os.makedirs(sim, exist_ok=True)
fdtd, csx, ports, nf, info, mesh = solver.build_simulation(p, sim)
if variant == "E2":
    x = np.asarray(mesh[0])
    k = int(np.argmin(np.abs(x - 87.4573)))
    x[k] = 87.4572999
    csx.GetGrid().SetLines("x", x)
x = np.asarray(csx.GetGrid().GetLines("x"))
print(variant, "x lines near port:", np.round(x[(x > 87.4) & (x < 87.6)], 7))
fdtd.Run(sim, cleanup=True, verbose=0)
i = np.loadtxt(os.path.join(sim, "port_it_1"), comments="%")
u = np.loadtxt(os.path.join(sim, "port_ut_1"), comments="%")
print(f"{variant}: max|u| {np.abs(u[:, 1]).max():.3e}  max|i| {np.abs(i[:, 1]).max():.3e}")
