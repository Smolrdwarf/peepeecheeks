"""Case 1 - straight WR-90 rectangular waveguide, TE10, 8-12 GHz.

HFSS equivalent: Driven Modal, two wave ports on an air-filled 22.86 x 10.16 mm box,
frequency sweep 8-12 GHz. HFSS reports S-parameters, port gamma (alpha + j*beta) and
port wave impedance; for an empty guide these converge to the analytic TE10 values.

openEMS: FDTD with two rectangular waveguide (mode-matched) ports and PML terminations.
beta is extracted numerically from the S21 phase difference of two guides that differ
only in length (removes any port reference-plane ambiguity).
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from common import CaseResult, Metric, sim_path, add_hfss_metrics, savefig

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0, Z0

CASE = "wr90_waveguide"
unit = 1e-3
a, b = 22.86, 10.16           # mm
f_start, f_stop = 8e9, 12e9
f_eval = 10e9
dL = 20.0                     # mm, length difference between the two runs


def run(length: float, tag: str):
    res = C0 / f_stop / unit / 30
    fdtd = openEMS(NrTS=30000, EndCriteria=1e-5)
    fdtd.SetGaussExcite(0.5 * (f_start + f_stop), 0.5 * (f_stop - f_start))
    fdtd.SetBoundaryCond([0, 0, 0, 0, 3, 3])
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(unit)
    mesh.AddLine('x', [0, a]); mesh.AddLine('y', [0, b]); mesh.AddLine('z', [0, length])
    p1 = ([0, 0, 10 * res], [a, b, 15 * res])
    p2 = ([0, 0, length - 10 * res], [a, b, length - 15 * res])
    for s, e in (p1, p2):
        mesh.AddLine('z', [s[2], e[2]])
    ports = [fdtd.AddRectWaveGuidePort(0, *p1, 'z', a * unit, b * unit, 'TE10', 1),
             fdtd.AddRectWaveGuidePort(1, *p2, 'z', a * unit, b * unit, 'TE10')]
    mesh.SmoothMeshLines('all', res, ratio=1.4)
    path = sim_path(f"{CASE}_{tag}")
    fdtd.Run(path, cleanup=True, verbose=0)
    freq = np.linspace(f_start, f_stop, 401)
    for p in ports:
        p.CalcPort(path, freq)
    s11 = ports[0].uf_ref / ports[0].uf_inc
    s21 = ports[1].uf_ref / ports[0].uf_inc
    zl = ports[0].uf_tot / ports[0].if_tot
    return freq, s11, s21, zl


def main() -> CaseResult:
    t0 = time.time()
    L1 = 100.0
    freq, s11, s21a, zl = run(L1, "L1")
    _, _, s21b, _ = run(L1 + dL, "L2")

    ph = np.angle(s21b / s21a)
    ph0 = np.mod(-ph[0], 2 * np.pi)                   # dphi(8 GHz) < 2*pi for dL = 20 mm
    dphi = ph0 + np.unwrap(-ph) - (-ph[0])
    beta_num = dphi / (dL * unit)

    k = 2 * np.pi * freq / C0
    kc = np.pi / (a * unit)
    beta_an = np.sqrt(k ** 2 - kc ** 2)
    zte_an = Z0 * k / beta_an
    fc_an = C0 / (2 * a * unit)
    # numerical cutoff from least-squares fit of beta^2 = k^2 - kc^2
    kc2_fit = np.mean(k ** 2 - beta_num ** 2)
    fc_num = np.sqrt(kc2_fit) * C0 / (2 * np.pi)

    i = np.argmin(abs(freq - f_eval))
    r = CaseResult(CASE, "WR-90 straight guide, TE10, 8-12 GHz",
                   "Driven Modal, 2 wave ports, interpolating sweep 8-12 GHz")
    r.metrics += [
        Metric("TE10 cutoff frequency", "GHz", fc_num / 1e9, fc_an / 1e9, "analytic", tol_pct=0.5),
        Metric("beta @ 10 GHz", "rad/m", beta_num[i], beta_an[i], "analytic", tol_pct=0.5),
        Metric("wave impedance @ 10 GHz", "ohm", float(np.real(zl[i])), zte_an[i], "analytic", tol_pct=1.0),
        Metric("min |S21| over band", "dB", float(20 * np.log10(abs(s21a).min())), 0.0, "analytic", tol_abs=0.1),
        Metric("max |S11| over band", "dB", float(20 * np.log10(abs(s11).max())), -30.0, "analytic", upper_bound=True),
    ]
    r.notes.append("max |S11|: an ideal matched guide has S11 -> -inf; treated as an upper limit of -30 dB.")
    add_hfss_metrics(r)
    r.runtime_s = time.time() - t0

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    ax[0].plot(freq / 1e9, 20 * np.log10(abs(s11)), label="S11"); ax[0].plot(freq / 1e9, 20 * np.log10(abs(s21a)), label="S21")
    ax[0].set(xlabel="f (GHz)", ylabel="dB", title="S-parameters (openEMS)"); ax[0].legend(); ax[0].grid()
    ax[1].plot(freq / 1e9, beta_an, "k-", label="analytic"); ax[1].plot(freq / 1e9, beta_num, "r--", label="openEMS")
    ax[1].set(xlabel="f (GHz)", ylabel="beta (rad/m)", title="Propagation constant"); ax[1].legend(); ax[1].grid()
    ax[2].plot(freq / 1e9, zte_an, "k-", label="analytic Z_TE"); ax[2].plot(freq / 1e9, np.real(zl), "r--", label="openEMS Re(Z)")
    ax[2].set(xlabel="f (GHz)", ylabel="ohm", title="Wave impedance"); ax[2].legend(); ax[2].grid()
    savefig(fig, CASE)
    return r


if __name__ == "__main__":
    r = main(); r.print(); r.save()
