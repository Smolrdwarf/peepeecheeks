"""Case 6 - the case-2 PEC cavity solved with a finite-element eigenmode solver (Elmer).

This is the open-source analogue of the HFSS Eigenmode solver: tetrahedral mesh (gmsh),
curl-conforming (Nedelec) edge elements of 2nd order (HFSS default order), generalized
eigenproblem  curl(1/mu curl E) = omega^2 eps E  solved directly - no time stepping and
no spectral peak picking as in the openEMS version.

Eigenvalues are omega^2 (SI). A shift-invert target of (2 pi 13 GHz)^2 keeps the
spurious lambda = 0 gradient modes further away than every physical mode below
~18.4 GHz, so the requested modes are all physical.
"""
import os, re, subprocess, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from common import CaseResult, Metric, sim_path, add_hfss_metrics
from cases.cavity_eigenmode import analytic_modes, a, b, d, f_lo, f_hi

CASE = "elmer_cavity_eigenmode"
ELMER_HOME = Path(os.environ.get("ELMER_HOME", Path.home() / "tools" / "elmer"))
H_MESH = 2.5e-3            # target tet edge length, m
F_SHIFT = 13e9
N_EIG = 12

SIF = """
Header
  Mesh DB "." "mesh"
End
Simulation
  Max Output Level = 5
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
  Active Solvers(1) = 1
End
Solver 1
  Equation = "VectorWave"
  Variable = E
  Procedure = "EMWaveSolver" "EMWaveSolver"
  Quadratic Approximation = Logical {quad}
  Linear System Solver = String "Direct"
  Linear System Direct Method = "umfpack"
  Use Global Mass Matrix = True
  Eigen Analysis = True
  Eigen System Values = {neig}
  Eigen System Convergence Tolerance = 1e-10
  Eigen System Select = smallest magnitude
  Eigen System Shift = Real {shift:.6e}
End
Boundary Condition 1
  Name = "PEC"
  Target Boundaries(1) = 1
  E {{e}} = Real 0
End
"""


def exe(name):
    p = ELMER_HOME / "bin" / f"{name}.exe"
    if not p.exists():
        raise FileNotFoundError(f"{p} not found - install Elmer (see README) or set ELMER_HOME")
    return str(p)


def make_mesh(work: Path, h: float):
    import gmsh
    gmsh.initialize(["", "-v", "0"])
    try:
        gmsh.model.add("cavity")
        vol = gmsh.model.occ.addBox(0, 0, 0, a * 1e-3, b * 1e-3, d * 1e-3)
        gmsh.model.occ.synchronize()
        surfs = [s[1] for s in gmsh.model.getBoundary([(3, vol)], oriented=False)]
        gmsh.model.addPhysicalGroup(2, surfs, 1, "PEC")
        gmsh.model.addPhysicalGroup(3, [vol], 1, "air")
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        gmsh.option.setNumber("Mesh.MeshSizeMin", h / 2)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.model.mesh.generate(3)
        n_tets = len(gmsh.model.mesh.getElementsByType(4)[0])
        gmsh.write(str(work / "mesh.msh"))
    finally:
        gmsh.finalize()
    subprocess.run([exe("ElmerGrid"), "14", "2", "mesh.msh", "-autoclean"], cwd=work,
                   check=True, capture_output=True, text=True)
    return n_tets


def solve(h=H_MESH, quadratic=True, tag=""):
    work = Path(sim_path(CASE + tag))
    for f in work.glob("*"):
        if f.is_file():
            f.unlink()
    n_tets = make_mesh(work, h)
    (work / "case.sif").write_text(SIF.format(quad=str(quadratic), neig=N_EIG, shift=(2 * np.pi * F_SHIFT) ** 2))
    (work / "ELMERSOLVER_STARTINFO").write_text("case.sif\n1\n")
    env = dict(os.environ, ELMER_HOME=str(ELMER_HOME), PATH=str(ELMER_HOME / "bin") + os.pathsep + os.environ["PATH"])
    p = subprocess.run([exe("ElmerSolver")], cwd=work, env=env, capture_output=True, text=True)
    (work / "solver.log").write_text(p.stdout + p.stderr)
    if p.returncode != 0:
        raise RuntimeError(f"ElmerSolver failed ({p.returncode}); see {work / 'solver.log'}")
    # EigenSolve prints lines like "EigenSolve:   1:  2.6826E+21  -1.2E+05"
    lam = [float(m.group(1)) for m in re.finditer(r"EigenSolve:\s+\d+:\s+([-+0-9.Ee]+)", p.stdout)]
    if not lam:
        raise RuntimeError(f"no eigenvalues parsed; see {work / 'solver.log'}")
    lam = np.array(lam)
    f = np.sqrt(lam[lam > 0]) / (2 * np.pi)
    return np.sort(f[f > 1e9]), n_tets


def main() -> CaseResult:
    t0 = time.time()
    freqs, n_tets = solve()
    r = CaseResult(CASE, f"PEC cavity {a} x {b} x {d} mm, FEM eigenmodes (Elmer, 2nd-order edge elements)",
                   "Eigenmode solver, PEC boundaries, N modes", solver="Elmer")
    for f_an, names in analytic_modes():
        j = np.argmin(abs(freqs - f_an)) if len(freqs) else None
        f_num = freqs[j] if j is not None and abs(freqs[j] - f_an) / f_an < 0.01 else np.nan
        r.metrics.append(Metric(f"f({'/'.join(names)})", "GHz", f_num / 1e9, f_an / 1e9, "analytic", tol_pct=0.5))
    r.notes.append(f"{n_tets:,} tetrahedra (h = {H_MESH*1e3:g} mm); eigenfrequencies found (GHz): "
                   + ", ".join(f"{x/1e9:.4f}" for x in freqs))
    add_hfss_metrics(r)
    r.runtime_s = time.time() - t0
    return r


if __name__ == "__main__":
    r = main(); r.print(); r.save()
