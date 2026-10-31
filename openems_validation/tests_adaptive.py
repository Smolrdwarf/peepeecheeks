"""Fast unit tests for adaptive.py (no solver): python tests_adaptive.py"""
import numpy as np
from adaptive import adaptive_mesh, max_delta_s

f = np.linspace(1, 2, 11)
# fake solver whose answer converges like 1/scale^2 toward S = exp(-j f)
fake = lambda s: {"S": np.exp(-1j * f) * (1 + 0.5 / s**2)}
quiet = lambda *_: None

r = adaptive_mesh(fake, start_scale=1, refine_factor=1.5, max_delta=0.02, log=quiet)
assert r.converged and r.passes[0].delta_s is None, r.table()
assert all(p.delta_s >= 0.02 for p in r.passes[1:-1]) and r.final.delta_s < 0.02

r2 = adaptive_mesh(fake, start_scale=1, refine_factor=1.5, max_delta=0.02, min_converged_passes=2, max_passes=10, log=quiet)
assert r2.converged and len(r2.passes) == len(r.passes) + 1

r3 = adaptive_mesh(fake, max_delta=1e-9, max_passes=3, log=quiet)
assert not r3.converged and len(r3.passes) == 3 and "max_passes" in r3.reason

r4 = adaptive_mesh(fake, estimate_cells=lambda s: int(1e5 * s**3), max_cells=300_000,
                   refine_factor=1.25, max_delta=1e-9, log=quiet)
assert not r4.converged and "max_cells" in r4.reason and all(p.cells <= 300_000 for p in r4.passes)

assert max_delta_s(np.array([1, 1j]), np.array([1, 0])) == 1.0
print("adaptive.py: all tests passed")
print(r.table())
