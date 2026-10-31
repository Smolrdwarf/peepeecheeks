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
- **Validation Check** with a mesh and memory estimate, **Analyze All** (F5) running in a subprocess
  with live progress and **Abort**
- **Results**: S-parameter plot with an m1 marker, Smith chart, input impedance,
  3D far-field directivity pattern, 2D E/H-plane cuts, Touchstone export
- Examples: half-wave dipole and probe-fed patch antenna (File > Examples)

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

The solver can also be run headless: `python -m fieldforge.solver project.ffproj outdir` (run from `hfss/`).

## Layout

| File | Purpose |
|---|---|
| `fieldforge/model.py` | Project data model, safe expression evaluator, material library, JSON save/load |
| `fieldforge/solver.py` | Project → CSXCAD/openEMS translation, auto-mesher, port/far-field post-processing |
| `fieldforge/viewer.py` | PyVista 3D modeler widget and picking |
| `fieldforge/plots.py` | Result viewers |
| `fieldforge/gui.py` | Main window, tree, properties, menus, solver process control |
| `fieldforge/examples.py` | Built-in example projects |

## Validation so far

| Case | Result | Reference |
|---|---|---|
| Half-wave dipole, 70 mm | D = 2.18 dBi, efficiency ≈ 100 % | 2.15 dBi theory |
| Patch 32×40 mm on RO4003C, 1.524 mm | S11 = −22 dB at 2.30 GHz, D = 6.6 dBi, efficiency 93 % | typical for a patch (6–7 dBi) |

## Known limitations / next steps

- FDTD is a staircased Cartesian mesh: curved surfaces are approximated, and thin features need a finer mesh
- Only lumped ports so far (no waveport); only one port excited per run, so you get one S-matrix column
- No boolean operations (unite/subtract), rotation or interactive mouse drawing yet
- No parametric sweep or optimizer yet (design variables are already in place for it)
- No field overlays (|E| on a plane) yet
