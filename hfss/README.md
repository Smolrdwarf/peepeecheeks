# FieldForge

A desktop 3D electromagnetic simulation app with an HFSS-style workflow:
draw geometry, assign materials and ports, set up an analysis, solve, and look at the results.
The solver is [openEMS](https://openems.de) (an FDTD method, whereas HFSS uses FEM). The GUI is
written from scratch in PySide6, PyVista/VTK and matplotlib. No Ansys code or files are used.

## Features

- **Project Manager tree**: Model (objects grouped by material), Boundaries, Excitations,
  Analysis, Results, Definitions (design variables, materials)
- **Properties panel**: every dimension is an expression that can use design variables (`-patch_W/2`)
- **3D Modeler**: box, cylinder, sphere, rectangle (sheet); click to select; region wireframe;
  undo/redo; duplicate/delete
- **Lumped ports**, **region boundaries** per face (PML / Mur / PEC / PMC), material library
  (FR4, Rogers, PTFE, alumina, …) and user materials
- **Adaptive passes** like HFSS (*Maximum Number of Passes*, *Maximum Delta S*, *Minimum Converged Passes*): the
  FDTD grid is rebuilt finer each pass until the S-parameters stop changing. Off by default; turned on for imported
  HFSS designs, using the HFSS project's own settings
- **Conductor-edge mesh refinement** (on by default), where FDTD loses accuracy first. See
  [Validation](#validation-so-far)
- **Validation Check** with a mesh and memory estimate, **Analyze All** (F5) running in a subprocess
  with live progress and **Abort**
- **Results**: S-parameter plot with an m1 marker, Smith chart, input impedance,
  3D far-field directivity pattern, 2D E/H-plane cuts, Touchstone export
- Examples: half-wave dipole and probe-fed patch antenna (File > Examples)

## Importing Ansys HFSS designs and comparing results

**File → Import Ansys HFSS Design (.aedt)** (Ctrl+I) reads `.aedt` and `.aedtz` projects. HFSS doesn't need to be
installed. If the project has several designs, you pick one. You get an **import report** listing exactly what was
translated, approximated or skipped.

| HFSS | FieldForge |
|---|---|
| Design and project variables (`$var`), with units | Design variables in model units. Still editable: changing one re-runs the HFSS modeler history |
| Modeler history: Box, Cylinder, Cone, Sphere, Rectangle, Circle, Ellipse, Polyline (line/arc/angular arc/spline, line/rectangle/circle cross-sections), CoverLines, Move, Rotate, Mirror, Unite, Subtract, Intersect, SweepAlongVector, SweepAroundAxis, SheetThicken, SweepFaceAlongNormal, BodyFromFace/Edge, Clone, DuplicateMirror/AlongLine/AroundAxis, Fillet, Chamfer, relative coordinate systems | Replayed with OpenCascade (`fieldforge/aedt/geometry.py`). Each final object becomes an imported (`cad`) object |
| Imported CAD bodies (`NativeBody`, stored as Parasolid) | Read from the original STEP file if it sits next to the `.aedt`. Otherwise faces/edges taken from them are rebuilt from the outlines HFSS stores (arcs recovered by circle fitting) |
| Materials | εr, μr, tan δ. σ > 1e5 S/m, or SolveInside = false, becomes PEC |
| Radiation boundary / region padding | Simulation region with PML on every face |
| Perfect E, Finite Conductivity | PEC (conductor loss not modelled) |
| Lumped Port | Lumped port spanning the port sheet along its integration line |
| Wave Port | **Approximation:** a lumped port along the integration line, narrowed to the signal conductor |
| Driven setup and sweep | Sweep range and point count. The solution frequency is used for the far field |
| Adaptive settings (`MaximumPasses`, `MaxDeltaS`, `MinimumConvergedPasses`) | Adaptive passes with the same limits. `PercentRefinement` has no FDTD equivalent: each pass re-meshes the whole grid 1.25× finer |

**Results → Load Reference Results** loads S-parameters exported from HFSS: Touchstone `.sNp`, or a report
exported as `.csv` (e.g. `dB(S(1,1))`). They're overlaid as dashed curves on the S-parameter plot and Smith chart.
**Results → Compare with Reference** shows both |S11| curves, their difference, and a metrics table: resonance
frequency, depth, −10 dB bandwidth, and RMS/max difference.

Tested against 70 real HFSS designs. All import without crashing. The only operations that don't reproduce are
some fillets/chamfers on edges whose IDs can't be matched; the corner is left sharp, and the import report says so.

Things to keep in mind when comparing:
- HFSS solves with adaptive FEM. FieldForge uses FDTD, starting at 20 cells/λ (adjustable in the setup); with
  adaptive passes it refines the whole grid uniformly instead of locally. Curved and diagonal edges are staircased.
- Conductor and surface losses are not modelled; dielectric loss is.
- FieldForge excites one port per run, so you get the first column of the S-matrix.

## Setup (Windows)

1. Download `openEMS_x64_v0.37.0_msvc.zip` from the
   [openEMS-Project releases](https://github.com/thliebig/openEMS-Project/releases) and unzip it to
   `~/tools/openEMS` (or set `OPENEMS_INSTALL_PATH` to wherever you put it).
2. Create a venv with a **python.org Python 3.13**. Anaconda's Python is not suitable, because its
   bundled Qt DLLs clash with PySide6.
   ```
   py -3.13 -m venv .venv
   .venv\Scripts\pip install %USERPROFILE%\tools\openEMS\python\csxcad-0.7.0-cp313-cp313-win_amd64.whl ^
                             %USERPROFILE%\tools\openEMS\python\openems-0.37.0-cp313-cp313-win_amd64.whl
   .venv\Scripts\pip install -r hfss\requirements.txt
   ```
3. Run `hfss\FieldForge.bat`, or `.venv\Scripts\python hfss\run.py [project.ffproj]`.

## Setup (macOS, Apple Silicon)

openEMS has no macOS release, so it is compiled from source. Tested on an M1 with macOS 26.

1. Install [Homebrew](https://brew.sh) and a **python.org Python 3.13**.
2. Run `hfss/setup_macos.sh`. It installs `cmake`, `boost` and `hdf5`, creates `.venv`, builds openEMS into
   `~/opt/openEMS` (or `OPENEMS_INSTALL_PATH`) and installs its Python bindings into `.venv`.
3. Double-click `hfss/FieldForge.command`, or run `.venv/bin/python hfss/run.py [project.ffproj]`.

The solver can also be run headless: `python -m fieldforge.solver project.ffproj outdir` (run from `hfss/`).

## Layout

| File | Purpose |
|---|---|
| `fieldforge/model.py` | Project data model, safe expression evaluator, material library, JSON save/load |
| `fieldforge/solver.py` | Project → CSXCAD/openEMS translation, auto-mesher, adaptive pass loop, port/far-field post-processing |
| `fieldforge/adaptive.py` | HFSS-style adaptive convergence loop (max delta-S between passes), solver-agnostic |
| `fieldforge/viewer.py` | PyVista 3D modeler widget and picking |
| `fieldforge/plots.py` | Result viewers |
| `fieldforge/gui.py` | Main window, tree, properties, menus, solver process control |
| `fieldforge/examples.py` | Built-in example projects |
| `fieldforge/aedt/parser.py` | Reader for the AEDT text project format (`.aedt`, `.aedtz`) |
| `fieldforge/aedt/units.py` | AEDT expression/unit conversion |
| `fieldforge/aedt/geometry.py` | OpenCascade replay of the HFSS modeler history, plus tessellation |
| `fieldforge/aedt/importer.py` | AEDT design → FieldForge project (materials, boundaries, ports, setup) and the import report |
| `fieldforge/cad.py` | Cached replay/tessellation for imported objects (rebuilt when variables change) |
| `fieldforge/reference.py` | Touchstone / HFSS CSV reader and comparison metrics |
| `tests/test_adaptive_import.py` | HFSS adaptive settings → FieldForge mapping, old project files still load (no solver needed) |
| `tests/microstrip_mesh_check.py` | eps_eff of a 50 Ω microstrip with each mesh option (needs openEMS, ~5 min) |
| `tests/patch_mesh_check.py` | Patch example resonance with edge refinement off/on/adaptive (needs openEMS, ~10 min) |

## Validation so far

| Case | Result | Reference |
|---|---|---|
| Half-wave dipole, 70 mm | D = 2.18 dBi, efficiency ≈ 100 % | 2.15 dBi theory |
| Patch 32×40 mm on RO4003C, 1.524 mm | S11 = −22 dB at 2.30 GHz, D = 6.6 dBi, efficiency 93 % | typical for a patch (6–7 dBi) |

macOS (M1, openEMS built from source at `d686f1d`):

| Case | Result | Solve time |
|---|---|---|
| Half-wave dipole, 70 mm | S11 = −24.2 dB at 1.855 GHz, D = 2.17 dBi, efficiency 101 % | about 30 s |
| Patch 32×40 mm on RO4003C, 1.524 mm | S11 = −22.9 dB at 2.30 GHz, D = 6.65 dBi, efficiency 87.7 % | about 100 s |

The patch efficiency differs from the Windows run above (93 %); the cause has not been established.

### Mesh options (validation-suite microstrip)

50 Ω microstrip, 1.1 mm wide on 0.508 mm RO4003C (εr 3.55 from the material library), lumped ports at both ends.
eps_eff comes from the S21 phase difference between a 30 mm and a 50 mm line. Reference: Hammerstad-Jensen, 2.779.
Source: `tests/microstrip_mesh_check.py`.

| Mesh | eps_eff error | Cells (50 mm line) |
|---|---|---|
| Fixed 20 cells/λ, no edge refinement (previous behaviour) | +7.0 % | 19k |
| Edge refinement (new default) | +0.95 % | 40k |
| Edge refinement + adaptive, converged in 2-3 passes | +0.6 % | 105k |

At 6 GHz the wavelength-based cell (1.36 mm) is wider than the 1.1 mm strip, so the edge cells are sized from the
conductor's own width (width/10) as well as the wavelength (cell/4).

### Solver fixes from this work

- **Sheet ports:** a lumped port with width (all imported HFSS lumped ports and approximated wave ports) now gets a
  mesh line through its centre. Without it, openEMS's voltage probe could fall between two lines
  ("only a 1D/line integration is allowed") and the run produced NaN S-parameters.
- **Mesh-line merging:** merging no longer averages two fixed planes together. Before, a thin substrate under a coarse
  grid lost both its ground plane (z=0) and its top surface; they were merged into a single line halfway between.
  The merge tolerance is now also capped at half the thin-layer and edge cell sizes.

## Known limitations / next steps

- FDTD is a staircased Cartesian mesh: curved surfaces are approximated, and thin features need a finer mesh.
  Adaptive passes refine globally, so each pass costs roughly 2× the previous one
- No convergence plot yet: adaptive pass history is saved in `results.npz` (`adaptive_history`) and printed to the log
- Only lumped ports natively. Imported wave ports are approximated. Only one port is excited per run, so you get one S-matrix column
- No boolean operations (unite/subtract), rotation or interactive mouse drawing yet
- No parametric sweep or optimizer yet (design variables are already in place for it)
- No field overlays (|E| on a plane) yet
