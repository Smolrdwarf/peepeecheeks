"""Case 7 - the case-1 WR-90 guide solved in the frequency domain with FEM (Elmer).

Open-source analogue of HFSS Driven Modal with two wave ports: VectorHelmholtz solver,
2nd-order Nedelec tetrahedra (gmsh), one solve per frequency. Ports are first-order
absorbing (Robin) conditions tuned to the TE10 propagation constant, with a TE10 surface
load on port 1 as the source - the same idea as an HFSS wave port, but with the port mode
given analytically instead of solved by a 2-D port eigenproblem.

S-parameters are extracted WITHOUT relying on the source normalization: the solved
E_y / sin(pi x / a) along the guide is least-squares fitted to a forward and a backward
wave with beta as a free parameter, so beta and S11 = backward/forward (at z=0) come from
the field solution itself. A wrong port BC would show up as reflection (large S11).

Time convention: Elmer uses exp(-i omega t), so the forward (+z) wave is exp(+j beta z).
Verified empirically: that term carries amplitude 1.006-1.027, matching the unit TE10
source on port 1; the exp(-j beta z) term is 70-80 dB smaller.
"""
import os, subprocess, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.optimize import minimize_scalar
from common import CaseResult, Metric, sim_path, add_hfss_metrics, savefig
from cases.elmer_cavity_eigenmode import exe, ELMER_HOME

CASE = "elmer_wr90_waveguide"
a, b, Lg = 22.86e-3, 10.16e-3, 100e-3    # m
H_MESH = 3e-3
C0 = 299792458.0
FREQS = np.linspace(8e9, 12e9, 9)
f_eval = 10e9

SIF = """
Header
  Mesh DB "." "mesh"
End
Simulation
  Max Output Level = 3
  Coordinate System = "Cartesian"
  Simulation Type = Steady State
  Steady State Max Iterations = 1
End
Constants
  Permittivity of Vacuum = 8.8541878128e-12
  Permeability of Vacuum = 1.25663706212e-6
End
Body 1
  Equation = 1
  Material = 1
End
Material 1
  Relative Permittivity = 1.0
  Relative Permeability = 1.0
End
Equation 1
  Active Solvers(3) = 1 2 3
  Angular Frequency = Real {w:.10e}
End
Solver 1
  Equation = "VectorHelmholtz"
  Procedure = "VectorHelmholtz" "VectorHelmholtzSolver"
  Quadratic Approximation = Logical True
  Use Piola Transform = Logical True
  Linear System Symmetric = False
  Linear System Solver = String "Direct"
  Linear System Direct Method = "umfpack"
End
Solver 2
  Equation = "calcfields"
  Procedure = "VectorHelmholtz" "VectorHelmholtzCalcFields"
  Calculate Elemental Fields = Logical False
  Calculate Nodal Fields = Logical True
  Calculate Electric field = Logical True
  Linear System Solver = "Direct"
  Linear System Direct Method = umfpack
End
Solver 3
  Equation = "result output"
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Output File Name = case
  Vtu Format = Logical True
  Ascii Output = Logical True
  Save Bulk Only = Logical True
End
! physical groups: 1 = port1 (z=0), 2 = port2 (z=L), 3 = PEC walls
Boundary Condition 1
  Target Boundaries(1) = 1
  Electric Robin Coefficient im = Real {beta:.10e}
  Magnetic Boundary Load Im 2 = Variable Coordinate 1
    Real MATC "2*{beta:.10e}*sin(pi*tx/{a})"
End
Boundary Condition 2
  Target Boundaries(1) = 2
  Electric Robin Coefficient im = Real {beta:.10e}
End
Boundary Condition 3
  Target Boundaries(1) = 3
  E re {{e}} = Real 0
  E im {{e}} = Real 0
End
"""


def make_mesh(work: Path, h: float):
    import gmsh
    gmsh.initialize(["", "-v", "0"])
    try:
        vol = gmsh.model.occ.addBox(0, 0, 0, a, b, Lg)
        gmsh.model.occ.synchronize()
        p1, p2, walls = [], [], []
        for dim, tag in gmsh.model.getBoundary([(3, vol)], oriented=False):
            z = gmsh.model.occ.getCenterOfMass(dim, tag)[2]
            (p1 if abs(z) < 1e-9 else p2 if abs(z - Lg) < 1e-9 else walls).append(tag)
        gmsh.model.addPhysicalGroup(2, p1, 1); gmsh.model.addPhysicalGroup(2, p2, 2)
        gmsh.model.addPhysicalGroup(2, walls, 3); gmsh.model.addPhysicalGroup(3, [vol], 1)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h); gmsh.option.setNumber("Mesh.MeshSizeMin", h / 2)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.model.mesh.generate(3)
        n = len(gmsh.model.mesh.getElementsByType(4)[0])
        gmsh.write(str(work / "mesh.msh"))
    finally:
        gmsh.finalize()
    subprocess.run([exe("ElmerGrid"), "14", "2", "mesh.msh", "-autoclean"], cwd=work, check=True, capture_output=True)
    return n


def solve_freq(work: Path, f: float):
    k0 = 2 * np.pi * f / C0
    beta = np.sqrt(k0 ** 2 - (np.pi / a) ** 2)
    (work / "case.sif").write_text(SIF.format(w=2 * np.pi * f, beta=beta, a=a))
    (work / "ELMERSOLVER_STARTINFO").write_text("case.sif\n1\n")
    for v in (work / "mesh").glob("case*.vtu"):
        v.unlink()
    env = dict(os.environ, ELMER_HOME=str(ELMER_HOME), PATH=str(ELMER_HOME / "bin") + os.pathsep + os.environ["PATH"])
    p = subprocess.run([exe("ElmerSolver")], cwd=work, env=env, capture_output=True, text=True)
    (work / f"solver_{f/1e9:.2f}GHz.log").write_text(p.stdout + p.stderr)
    if p.returncode != 0:
        raise RuntimeError(f"ElmerSolver failed at {f/1e9} GHz; see {work}")
    return read_vtu(next((work / "mesh").glob("case*.vtu")))


def read_vtu(path: Path):
    import meshio
    m = meshio.read(path)
    return m.points, m.point_data


def fit_waves(pts, pdata, beta_guess):
    """Fit E_y / sin(pi x/a) = F e^{+j beta z} + R e^{-j beta z} (Elmer: exp(-i w t), so F is the
    +z-travelling wave). Returns beta, F (forward), R (reflected), rel. residual, cross-pol ratio."""
    E = pdata["electric field re"] + 1j * pdata["electric field im"]
    sx = np.sin(np.pi * pts[:, 0] / a)
    sel = sx > 0.3
    z, u = pts[sel, 2], E[sel, 1] / sx[sel]

    def lsq(beta):
        M = np.c_[np.exp(1j * beta * z), np.exp(-1j * beta * z)]
        c, *_ = np.linalg.lstsq(M, u, rcond=None)
        return c, np.linalg.norm(M @ c - u) / np.linalg.norm(u)

    grid = np.linspace(0.6 * beta_guess, 1.4 * beta_guess, 401)     # coarse scan, then refine
    b0 = grid[np.argmin([lsq(bb)[1] for bb in grid])]
    opt = minimize_scalar(lambda bb: lsq(bb)[1], bounds=(b0 - 0.01 * beta_guess, b0 + 0.01 * beta_guess),
                          method="bounded", options={"xatol": 1e-6})
    (A, B), res = lsq(opt.x)
    purity = np.linalg.norm(E[:, [0, 2]]) / np.linalg.norm(E)        # |Ex,Ez| / |E|: 0 for pure TE10
    return opt.x, A, B, res, purity


def main() -> CaseResult:
    t0 = time.time()
    work = Path(sim_path(CASE))
    n_tets = make_mesh(work, H_MESH)
    k0 = 2 * np.pi * FREQS / C0
    beta_an = np.sqrt(k0 ** 2 - (np.pi / a) ** 2)
    beta_num, s11, resid, purity, amp = [], [], [], [], []
    for f, bg in zip(FREQS, beta_an):
        pts, pdata = solve_freq(work, f)
        bt, A, B, res, pur = fit_waves(pts, pdata, bg)
        beta_num.append(bt); s11.append(B / A); resid.append(res); purity.append(pur); amp.append(abs(A))
    beta_num, s11 = np.array(beta_num), np.array(s11)
    fc_num = np.sqrt(np.mean(k0 ** 2 - beta_num ** 2)) * C0 / (2 * np.pi)
    i = np.argmin(abs(FREQS - f_eval))

    r = CaseResult(CASE, "WR-90 straight guide, TE10, 8-12 GHz, FEM frequency domain (Elmer)",
                   "Driven Modal, 2 wave ports, discrete sweep 8-12 GHz", solver="Elmer")
    r.metrics += [
        Metric("TE10 cutoff frequency", "GHz", fc_num / 1e9, C0 / (2 * a) / 1e9, "analytic", tol_pct=0.5),
        Metric("beta @ 10 GHz", "rad/m", beta_num[i], beta_an[i], "analytic", tol_pct=0.5),
        Metric("max |S11| over band", "dB", 20 * np.log10(np.abs(s11).max()), -30.0, "analytic", upper_bound=True),
    ]
    r.notes.append(f"{n_tets:,} tetrahedra (h = {H_MESH*1e3:g} mm), 2nd-order edge elements, {len(FREQS)} frequencies")
    r.notes.append(f"fit residual max {max(resid):.2e}; cross-polar |Ex,Ez|/|E| max {max(purity):.2e}; "
                   f"incident amplitude |A| {min(amp):.3f}-{max(amp):.3f} (1.0 expected from the port load)")
    add_hfss_metrics(r)
    r.runtime_s = time.time() - t0

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(FREQS / 1e9, beta_an, "k-", label="analytic"); ax[0].plot(FREQS / 1e9, beta_num, "bo", label="Elmer FEM")
    ax[0].set(xlabel="f (GHz)", ylabel="beta (rad/m)", title="Propagation constant"); ax[0].legend(); ax[0].grid()
    ax[1].plot(FREQS / 1e9, 20 * np.log10(np.abs(s11)), "bo-", label="S11 (Elmer)")
    ax[1].set(xlabel="f (GHz)", ylabel="dB", title="Reflection (port absorption quality)"); ax[1].legend(); ax[1].grid()
    savefig(fig, CASE)
    return r


if __name__ == "__main__":
    r = main(); r.print(); r.save()
