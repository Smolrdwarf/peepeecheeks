# Third-party licences

FieldForge is licensed under GPL-3.0-or-later (see `LICENSE`). It uses the components below. None of them is vendored in this repository. They are installed separately, and each keeps its own licence.

GPL-3.0-or-later is required in practice because FieldForge imports openEMS, which is GPL-3.0-or-later. Every other licence below is compatible with it.

| Component | Used for | Licence | Checked against |
|---|---|---|---|
| [openEMS](https://github.com/thliebig/openEMS) | FDTD solver | GPL-3.0-or-later | Installed package metadata (0.37.0) |
| [CSXCAD](https://github.com/thliebig/CSXCAD) | Geometry and materials for openEMS | LGPL-3.0-or-later | Installed package metadata (0.7.0) |
| [Elmer FEM](https://github.com/ElmerCSC/elmerfem) | FEM solver, run as a separate program | GPL-2.0, with the ElmerSolver library under LGPL-2.1 | Project repository |
| [Gmsh](https://gmsh.info) | Tetrahedral meshing for the Elmer cases | GPL-2.0-or-later | PyPI metadata (4.15.2) |
| [PySide6](https://doc.qt.io/qtforpython/) | Desktop GUI | LGPL-3.0-only OR GPL-2.0-only OR GPL-3.0-only | Installed package metadata (6.12.0) |
| [cadquery-ocp](https://github.com/CadQuery/OCP) | Python bindings for OpenCascade, used to replay `.aedt` geometry | Apache-2.0 | PyPI metadata (8.0.1.1.0) |
| [Open CASCADE Technology](https://dev.opencascade.org/resources/licensing) | CAD kernel bundled in cadquery-ocp | LGPL-2.1 with Open CASCADE exception | Open CASCADE licensing page |
| [PyVista](https://github.com/pyvista/pyvista) | 3D viewport | MIT | Installed package metadata (0.49.2) |
| [pyvistaqt](https://github.com/pyvista/pyvistaqt) | PyVista and Qt integration | MIT | Installed package metadata (0.13.1) |
| [Matplotlib](https://matplotlib.org) | Plots | Matplotlib License (PSF-based) | Installed package metadata (3.11.2) |
| [NumPy](https://numpy.org) | Arrays | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0 | Installed package metadata (2.5.4) |
| [SciPy](https://scipy.org) | Validation analysis | BSD-3-Clause | PyPI metadata (1.18.1) |
| [h5py](https://www.h5py.org) | Reading openEMS output | BSD-3-Clause | Installed package metadata (3.16.0) |
| [meshio](https://github.com/nschloe/meshio) | Reading Elmer output | MIT | PyPI metadata (5.3.5) |

Ansys, HFSS and AEDT are trademarks of their owners. This project contains no Ansys code and is not affiliated with Ansys or Synopsys.
