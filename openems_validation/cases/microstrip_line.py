"""Case 3 - 50 ohm microstrip on Rogers RO4003C (er = 3.38, h = 0.508 mm).

HFSS equivalent: Driven Modal/Terminal, wave ports at both ends of a straight line,
reports port Zo and gamma (eps_eff = (beta/k0)^2). Reference here: Hammerstad-Jensen
closed-form (accurate to ~0.2-1% for this geometry), evaluated at 2 GHz where
dispersion is negligible for h = 0.508 mm.

openEMS: MSL ports (3 voltage / 2 current probes) give numerical Z0 and beta.

Mesh convergence (2026-10-31, openEMS v0.37): the openEMS-tutorial mesh (lambda/50,
edge cells res/4, 5 substrate cells) gave eps_eff +8.6% / Z0 -2.6%. Refining the strip
edges (res/8) and substrate (10 cells) at lambda/80 gives +0.5% / -0.6%. Defaults below
are the converged settings; eps_eff is dominated by the strip-edge cell size.
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from common import CaseResult, Metric, sim_path, add_hfss_metrics, savefig

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0, Z0 as ETA0

CASE = "microstrip_line"
unit = 1e-6                  # um
er, h, w = 3.38, 508.0, 1100.0
L = 30000.0                  # half-length of the line
f_max, f_eval = 6e9, 2e9


def hammerstad_jensen(er, w, h):
    u = w / h
    a = 1 + np.log((u**4 + (u / 52) ** 2) / (u**4 + 0.432)) / 49 + np.log(1 + (u / 18.1) ** 3) / 18.7
    b = 0.564 * ((er - 0.9) / (er + 3)) ** 0.053
    eeff = (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / u) ** (-a * b)
    fu = 6 + (2 * np.pi - 6) * np.exp(-((30.666 / u) ** 0.7528))
    z01 = ETA0 / (2 * np.pi) * np.log(fu / u + np.sqrt(1 + (2 / u) ** 2))
    return z01 / np.sqrt(eeff), eeff


def main(res_div=80, edge_div=8, nz_sub=11, top=12, side=15, tag="") -> CaseResult:
    t0 = time.time()
    fdtd = openEMS(EndCriteria=1e-5)
    fdtd.SetGaussExcite(f_max / 2, f_max / 2)
    fdtd.SetBoundaryCond(['PML_8', 'PML_8', 'MUR', 'MUR', 'PEC', 'MUR'])
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(unit)

    res = C0 / (f_max * np.sqrt(er)) / unit / res_div
    third = np.array([2 * res / 3, -res / 3]) / edge_div          # "thirds rule" at strip edges
    mesh.AddLine('x', [-L, 0, L]); mesh.SmoothMeshLines('x', res)
    mesh.AddLine('y', 0); mesh.AddLine('y', w / 2 + third); mesh.AddLine('y', -w / 2 - third)
    mesh.SmoothMeshLines('y', res / edge_div)
    mesh.AddLine('y', [-side * w, side * w]); mesh.SmoothMeshLines('y', res)
    mesh.AddLine('z', np.linspace(0, h, nz_sub)); mesh.AddLine('z', top * h); mesh.SmoothMeshLines('z', res)

    sub = csx.AddMaterial('RO4003C', epsilon=er)
    sub.AddBox([-L, -side * w, 0], [L, side * w, h])
    pec = csx.AddMetal('PEC')
    ports = [
        fdtd.AddMSLPort(1, pec, [-L, -w / 2, h], [0, w / 2, 0], 'x', 'z', excite=-1,
                        FeedShift=10 * res, MeasPlaneShift=L / 3, priority=10),
        fdtd.AddMSLPort(2, pec, [L, -w / 2, h], [0, w / 2, 0], 'x', 'z',
                        MeasPlaneShift=L / 3, priority=10),
    ]
    path = sim_path(CASE + tag)
    fdtd.Run(path, cleanup=True, verbose=0)

    freq = np.linspace(0.5e9, f_max, 221)
    ports[0].CalcPort(path, freq)                 # Z_ref = numerical line impedance from probes
    z0_num = np.real(ports[0].Z_ref)
    eeff_num = (np.real(ports[0].beta) / (2 * np.pi * freq / C0)) ** 2
    for p in ports:
        p.CalcPort(path, freq, ref_impedance=50)  # S-parameters referenced to 50 ohm
    s11 = ports[0].uf_ref / ports[0].uf_inc
    s21 = ports[1].uf_ref / ports[0].uf_inc
    z0_an, eeff_an = hammerstad_jensen(er, w, h)
    np.savetxt(Path(path) / "curves.csv", np.c_[freq, z0_num, eeff_num], delimiter=",", header="f,z0,eeff")
    i = np.argmin(abs(freq - f_eval))

    r = CaseResult(CASE, f"Microstrip w={w/1e3} mm on h={h/1e3} mm RO4003C (er={er})",
                   "Driven Terminal, 2 wave ports, sweep 0.5-6 GHz (port Zo, gamma)")
    r.metrics += [
        Metric("Z0 @ 2 GHz", "ohm", z0_num[i], z0_an, "analytic", tol_pct=2.0),
        Metric("eps_eff @ 2 GHz", "-", eeff_num[i], eeff_an, "analytic", tol_pct=2.0),
        Metric("min |S21| 0.5-6 GHz (lossless)", "dB", 20 * np.log10(abs(s21).min()), 0.0, "analytic", tol_abs=0.2),
    ]
    r.notes.append("Reference is Hammerstad-Jensen (zero strip thickness, quasi-static).")
    add_hfss_metrics(r)
    r.runtime_s = time.time() - t0

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    ax[0].plot(freq / 1e9, z0_num, "r-", label="openEMS"); ax[0].axhline(z0_an, color="k", ls="--", label="Hammerstad-Jensen")
    ax[0].set(xlabel="f (GHz)", ylabel="ohm", title="Characteristic impedance"); ax[0].legend(); ax[0].grid()
    ax[1].plot(freq / 1e9, eeff_num, "r-", label="openEMS"); ax[1].axhline(eeff_an, color="k", ls="--", label="Hammerstad-Jensen")
    ax[1].set(xlabel="f (GHz)", ylabel="eps_eff", title="Effective permittivity"); ax[1].legend(); ax[1].grid()
    ax[2].plot(freq / 1e9, 20 * np.log10(abs(s11)), label="S11 (50 ohm ref)"); ax[2].plot(freq / 1e9, 20 * np.log10(abs(s21)), label="S21")
    ax[2].set(xlabel="f (GHz)", ylabel="dB", title="S-parameters"); ax[2].legend(); ax[2].grid()
    savefig(fig, CASE)
    return r


if __name__ == "__main__":
    r = main(); r.print(); r.save()
