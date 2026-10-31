"""End-to-end check of FieldForge's mesh options on the validation-suite microstrip (case 3).

50 ohm microstrip (w = 1.1 mm on 0.508 mm RO4003C), lumped ports at both ends. eps_eff is taken
from the S21 phase difference between a 30 mm and a 50 mm line (port parasitics cancel), and
compared for: the original mesh, conductor-edge refinement, and adaptive passes.
Reference: Hammerstad-Jensen with the library's eps_r for RO4003C (3.55 -> eps_eff 2.779).

Run:  python hfss/tests/microstrip_mesh_check.py    (~5-10 min, needs openEMS)
"""
import os
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fieldforge.model import LumpedPort, Primitive, Project, Region, Setup   # noqa: E402
from fieldforge.solver import solve                                          # noqa: E402

C0 = 299792458.0
W, H, L_SHORT, L_LONG, F_EVAL = 1.1, 0.508, 30.0, 50.0, 2e9


def hammerstad_jensen_eeff(er, w, h):
    u = w / h
    a = 1 + np.log((u**4 + (u / 52) ** 2) / (u**4 + 0.432)) / 49 + np.log(1 + (u / 18.1) ** 3) / 18.7
    b = 0.564 * ((er - 0.9) / (er + 3)) ** 0.053
    return (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / u) ** (-a * b)


ER = Project().material("Rogers RO4003C").epsr            # FieldForge library value (3.55)
EEFF_HJ = hammerstad_jensen_eeff(ER, W, H)


def microstrip(length, **setup_kw):
    p = Project(name=f"MSL{length:g}", units="mm")
    p.variables = {"w": str(W), "h": str(H), "L": str(length)}
    side, margin = 15 * W, 5.0   # ports must not sit on the absorbing region boundary
    p.objects = [
        Primitive("Substrate", "box", "Rogers RO4003C",
                  params={"x": f"-L/2-{margin}", "y": f"-{side}", "z": "0", "dx": f"L+{2 * margin}",
                          "dy": f"{2 * side}", "dz": "h"}),
        Primitive("Strip", "rectangle", "pec", axis="Z",
                  params={"x": "-L/2", "y": "-w/2", "z": "h", "w": "L", "l": "w"}),
    ]
    p.ports = [LumpedPort("P1", axis="Z", params={"x1": "-L/2", "y1": "-w/2", "z1": "0",
                                                 "x2": "-L/2", "y2": "w/2", "z2": "h"}),
               LumpedPort("P2", axis="Z", excited=False,
                          params={"x1": "L/2", "y1": "-w/2", "z1": "0", "x2": "L/2", "y2": "w/2", "z2": "h"})]
    # ground = PEC zmin face; open (MUR) elsewhere
    p.region = Region(boundaries=["MUR", "MUR", "MUR", "MUR", "PEC", "MUR"],
                      explicit=[-length / 2 - margin, length / 2 + margin, -side, side, 0.0, 12 * H])
    p.setup = Setup(f_start="0.5", f_stop="6", n_points=221, farfield=False, end_criteria_db=-50, **setup_kw)
    return p


def run(project, work):
    pf = os.path.join(work, project.name + ".ffproj")
    project.save(pf)
    d = np.load(solve(pf, os.path.join(work, project.name)))
    return d


def eeff(cfg, work):
    out = []
    for L in (L_SHORT, L_LONG):
        out.append(run(microstrip(L, **cfg), work))
    f = out[0]["freq"]
    ratio = out[1]["S"][1] / out[0]["S"][1]                       # S21(long) / S21(short)
    dphi = np.unwrap(-np.angle(ratio))
    dphi += 2 * np.pi * np.round((0 - dphi[0]) / (2 * np.pi))      # ~0 rad at the lowest frequency
    beta = dphi / ((L_LONG - L_SHORT) * 1e-3)
    e = (beta / (2 * np.pi * f / C0)) ** 2
    hist = [o["adaptive_history"] for o in out if "adaptive_history" in o]
    return float(np.interp(F_EVAL, f, e)), [len(o["mesh_x"]) * len(o["mesh_y"]) * len(o["mesh_z"]) for o in out], hist


if __name__ == "__main__":
    work = os.path.join(tempfile.gettempdir(), "fieldforge_msl_check")
    os.makedirs(work, exist_ok=True)
    configs = {
        "original mesh (no edge refinement)": dict(edge_refinement=False),
        "edge refinement": dict(edge_refinement=True),
        "edge refinement + adaptive (dS < 0.02)": dict(edge_refinement=True, adaptive=True),
    }
    rows = []
    for name, cfg in configs.items():
        t0 = time.time()
        e, cells, hist = eeff(cfg, work)
        rows.append((name, e, cells, time.time() - t0, hist))
        print(f"\n>>> {name}: eps_eff = {e:.4f} ({100 * (e / EEFF_HJ - 1):+.2f} % vs HJ), cells {cells}", flush=True)
    print("\n| Mesh | eps_eff @ 2 GHz | vs Hammerstad-Jensen | cells (30 / 50 mm) | time |")
    print("|---|---|---|---|---|")
    for name, e, cells, t, hist in rows:
        print(f"| {name} | {e:.4f} | {100 * (e / EEFF_HJ - 1):+.2f} % | {cells[0]:,} / {cells[1]:,} | {t:.0f} s |")
        for h in hist:
            print("|  passes: " + "; ".join(f"#{int(r[0])} x{r[1]:.2f} {int(r[2]):,} cells dS={r[3]:.3f}" for r in h) + " | | | | |")
