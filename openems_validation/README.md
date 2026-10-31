# openEMS validation suite

Benchmark cases that check the open-source solvers (openEMS = FDTD, Elmer = FEM) against the results
Ansys HFSS would give for the same structures. Each case has a defined "HFSS equivalent"
setup so the numbers can be compared one-to-one.

**Status (2026-10-31):** Ansys HFSS is not installed on the dev machine, so every case is
currently compared against an **analytic / closed-form reference** (the same values HFSS
converges to for these canonical structures). No HFSS numbers have been compared yet. To add
them, see [HFSS reference data](#hfss-reference-data).

Latest results: [`results/SUMMARY.md`](results/SUMMARY.md)

| Case | HFSS equivalent | Reference | Result |
|---|---|---|---|
| `wr90_waveguide` | Driven Modal, 2 wave ports, 8-12 GHz | TE10 analytic (fc, beta, Z_TE) | PASS (≤0.6 %) |
| `cavity_eigenmode` | Eigenmode solver, PEC box | analytic f_mnp, 8 modes | PASS (≤0.05 %) |
| `microstrip_line` | Driven Terminal, port Zo / gamma | Hammerstad-Jensen | PASS (Z0 −0.6 %, eps_eff +0.5 %) |
| `microstrip_adaptive` | Adaptive passes, Maximum Delta S 0.02 | Hammerstad-Jensen + delta-S criterion | PASS (converges in 3 passes from the tutorial mesh) |
| `elmer_cavity_eigenmode` | Eigenmode solver (FEM, like HFSS) | analytic f_mnp, 8 modes | PASS (≤0.009 %) |
| `elmer_wr90_waveguide` | Driven Modal, wave ports (FEM) | TE10 analytic | PASS (beta −0.014 %, S11 −67 dB) |
| `pec_sphere_rcs` | Plane-wave incidence, monostatic RCS | Mie series | **FAIL** at ka=3 (−0.79 dB vs ±0.5 dB); RMS 0.40 dB passes |

## openEMS (FDTD) vs Elmer (FEM)

HFSS is an FEM code, so Elmer is the closer analogue. Elmer runs here with gmsh tetrahedra and second-order
Nedelec edge elements, which is HFSS's default order. On the same structures:

| | openEMS (FDTD) | Elmer (FEM) |
|---|---|---|
| Cavity, worst mode error | 0.05 % | 0.009 % |
| Cavity runtime | 44 s (150k timesteps, peak picking) | 1.7 s (direct eigen solve) |
| Cavity degenerate TE111/TM111 | one peak | two modes, resolved |
| WR-90 beta @ 10 GHz | +0.07 % | −0.014 % |
| WR-90 max S11 (port quality) | −47.7 dB | −66.8 dB |
| WR-90 runtime | 9 s for the whole band | 62 s (9 frequencies × ~7 s) |

Takeaways:
- FEM is the better fit for eigenmodes and for narrowband or high-accuracy work.
- FDTD is cheaper for broadband sweeps and large open problems.
- A product that mirrors HFSS would plausibly use FEM by default (Driven Modal and Eigenmode), with FDTD
  available for transient analysis.

Elmer caveats:
- **Port mode is analytic.** The wave-port boundary uses the analytic TE10 mode and beta. HFSS solves a 2D port
  eigenproblem instead; Elmer can do that too (see its `EM_port_eigen` tests), but it isn't wired up yet.
- **Different time convention.** Elmer uses exp(−iωt), the opposite of the usual engineering convention, so
  phases flip sign relative to HFSS/Touchstone.
- **Palace not included.** AWS's FEM solver Palace was the other candidate. It has no Windows build: it needs
  WSL plus a Spack compile, so it was deferred.

## Install (Windows)

```powershell
powershell -ExecutionPolicy Bypass -File install_openems.ps1      # -> ~/tools/openEMS + ~/tools/openems-venv
```

The script downloads the openEMS v0.37.0 MSVC build and Elmer 26.1 (native Windows, no MPI; set
`ELMER_HOME` to use another location). It unpacks both, creates a venv, and installs gmsh/meshio and
the bundled CSXCAD/openEMS wheels (they are not on PyPI; Python 3.13 or 3.14 required).
`common.py` points `CSXCAD_INSTALL_PATH` at `~/tools/openEMS` unless it is already set.

## Run

```bash
~/tools/openems-venv/Scripts/python run_all.py                  # all cases (~5 min on 16 cores)
~/tools/openems-venv/Scripts/python run_all.py microstrip_line  # one case
```

Each case writes `results/<case>.json` and `results/<case>.png`; `run_all.py` also writes
`results/SUMMARY.md`. Raw solver output goes to `%TEMP%/openems_validation/<case>` (not in git).
Peak memory is ~1-2 GB (sphere case, 2.9 M cells).

## Findings so far

- **The openEMS tutorial mesh isn't good enough on microstrip.** It gave eps_eff +8.6 %. You need strip-edge cells of
  about res/8 and 10 cells across the substrate to get within 0.5 %. HFSS hides this with adaptive meshing; our
  tool will need an equivalent convergence loop.
- **Sphere RCS has an error from the stepped (staircased) surface.** A 12 mm uniform grid gives 0.83 dB RMS; a
  6 mm region around the sphere gives 0.40 dB. The remaining error sits next to the deep nulls. A 4 mm region
  (6.1 M cells, 6.5 min) gives 0.34 dB RMS but leaves ka=3 unchanged at −0.77 dB. So the ka=3 error doesn't
  come from how finely the sphere is meshed. Unverified suspects: grid dispersion in the 12 mm background cells,
  or the plane-wave/NF2FF box setup. 6 mm stays the default.
- **Watch for stale far-field results.** `nf2ff.CalcNF2FF` reuses an existing result file by default
  (`read_cached=True`). That silently returned stale results after a mesh change; the sphere case passes
  `read_cached=False`.
- **openEMS has no eigenmode solver.** The cavity case approximates one by recording a long ringing time signal
  and finding the spectral peaks. That works for lossless cavities but won't give Q factors the way HFSS does.

## Adaptive mesh convergence (`adaptive.py`)

The loop now lives in FieldForge (`hfss/fieldforge/adaptive.py`), which uses it for its adaptive passes;
`adaptive.py` here re-exports it, so this folder needs the `hfss/` folder next to it.

HFSS refines its mesh until S-parameters change by less than *Maximum Delta S* between passes.
openEMS has no adaptive meshing, so `adaptive.py` does the same from the outside:

```python
from adaptive import adaptive_mesh
res = adaptive_mesh(simulate,            # scale -> {"S": complex array (n_freq, ...), ...}
                    estimate_cells,      # scale -> int, checked against max_cells before each pass
                    start_scale=1.0, refine_factor=1.25, max_delta=0.02,
                    min_converged_passes=1, max_passes=6, max_cells=5_000_000)
print(res.table()); best = res.final.output
```

- Delta S is the max of |S_n − S_(n−1)| over **all frequencies** and all S-matrix entries. HFSS uses a
  single adaptive frequency, so this is stricter.
- Refinement is global (one mesh-density scale), not error-driven local refinement like HFSS.
- Each pass is a full rerun. Cost grows roughly with scale^3 to scale^4.
- Like HFSS, it can stop too early if two coarse passes happen to agree.
  `min_converged_passes=2` is the safer setting for unfamiliar geometry.

On the microstrip, it starts from the tutorial mesh (eps_eff +8.6 %). Pass 2 changes S by
0.25, pass 3 by 0.0014, and it stops at Z0 −0.6 % / eps_eff +0.5 %. Unit tests that don't need the solver:
`python tests_adaptive.py`.

## HFSS reference data

Put HFSS results in `references/hfss/<case>.csv` as `metric,value` rows. Metric names must
match the `name` column in `results/<case>.json`, in the units listed there. Example:

```
metric,value
TE10 cutoff frequency,6.5571
beta @ 10 GHz,158.24
```

When the file exists, each matching metric is compared against HFSS as well as the analytic value,
using the same tolerance. Use the geometry and settings documented at the top of each `cases/*.py`.

## Layout

```
adaptive.py        re-exports hfss/fieldforge/adaptive.py (single implementation, also used by FieldForge)
tests_adaptive.py  solver-free unit tests for adaptive.py
common.py          Metric/CaseResult, tolerance checks, HFSS CSV loader, paths
run_all.py         runs cases, writes results/SUMMARY.md, exit code 1 on any failure
cases/*.py         one benchmark per file; each runnable on its own (elmer_* use Elmer + gmsh)
references/hfss/   HFSS exports (empty for now)
results/           JSON + plots + SUMMARY.md from the last run
```
