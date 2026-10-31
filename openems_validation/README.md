# openEMS validation suite

Benchmark cases that check the open-source solver (openEMS, FDTD) against the results
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
| `pec_sphere_rcs` | Plane-wave incidence, monostatic RCS | Mie series | **FAIL** at ka=3 (−0.79 dB vs ±0.5 dB); RMS 0.40 dB passes |

## Install (Windows)

```powershell
powershell -ExecutionPolicy Bypass -File install_openems.ps1      # -> ~/tools/openEMS + ~/tools/openems-venv
```

The script downloads the openEMS v0.37.0 MSVC build, unpacks it, creates a venv and installs
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
  6 mm region around the sphere gives 0.40 dB. The remaining error sits next to the deep nulls. A 4 mm run was
  stopped because the machine was low on memory, so convergence below 6 mm hasn't been verified.
- **Watch for stale far-field results.** `nf2ff.CalcNF2FF` reuses an existing result file by default
  (`read_cached=True`). That silently returned stale results after a mesh change; the sphere case passes
  `read_cached=False`.
- **openEMS has no eigenmode solver.** The cavity case approximates one by recording a long ringing time signal
  and finding the spectral peaks. That works for lossless cavities but won't give Q factors the way HFSS does.

## Adaptive mesh convergence (`adaptive.py`)

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
adaptive.py        HFSS-style adaptive mesh loop (max delta-S convergence)
tests_adaptive.py  solver-free unit tests for adaptive.py
common.py          Metric/CaseResult, tolerance checks, HFSS CSV loader, paths
run_all.py         runs cases, writes results/SUMMARY.md, exit code 1 on any failure
cases/*.py         one benchmark per file; each runnable on its own
references/hfss/   HFSS exports (empty for now)
results/           JSON + plots + SUMMARY.md from the last run
```
