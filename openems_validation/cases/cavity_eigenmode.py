"""Case 2 - closed rectangular PEC cavity, resonant frequencies.

HFSS equivalent: Eigenmode solution type, air box with PEC walls, request N modes.
HFSS returns complex eigenfrequencies (purely real for a lossless PEC cavity).

openEMS has no eigenmode solver, so we ring the cavity: broadband soft E-field kick at an
asymmetric point, record E-field voltage probes for a long time, then pick spectral peaks
(Hann window + zero-padded FFT + parabolic peak interpolation).
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import itertools
import numpy as np
from common import CaseResult, Metric, sim_path, add_hfss_metrics, savefig

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0

CASE = "cavity_eigenmode"
unit = 1e-3
a, b, d = 22.86, 10.16, 30.0       # mm (x, y, z)
f_lo, f_hi = 5e9, 17.5e9
cell = 0.5                          # mm
NrTS = 150000


def analytic_modes():
    modes = {}
    for m, n, p in itertools.product(range(5), range(3), range(5)):
        te = p >= 1 and (m, n) != (0, 0)
        tm = m >= 1 and n >= 1
        if not (te or tm):
            continue
        f = C0 / 2 * np.sqrt((m / (a * unit)) ** 2 + (n / (b * unit)) ** 2 + (p / (d * unit)) ** 2)
        if f_lo <= f <= f_hi:
            names = ([f"TE{m}{n}{p}"] if te else []) + ([f"TM{m}{n}{p}"] if tm else [])
            key = round(f / 1e6, 3)
            modes.setdefault(key, (f, []))[1].extend(names)
    return sorted(modes.values())


def snap(lines, v):
    return lines[np.argmin(abs(lines - v))]


def main() -> CaseResult:
    t0 = time.time()
    fdtd = openEMS(NrTS=NrTS, EndCriteria=0)
    fdtd.SetGaussExcite(0.5 * (f_lo + 20e9), 0.5 * (20e9 - f_lo))
    fdtd.SetBoundaryCond([0, 0, 0, 0, 0, 0])
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(unit)
    L = {}
    for ax, size in zip("xyz", (a, b, d)):
        L[ax] = np.linspace(0, size, int(np.ceil(size / cell)) + 1)
        mesh.SetLines(ax, L[ax])

    # excitation: one-cell box at an asymmetric spot, all three E components
    p0 = [snap(L[ax], v) for ax, v in zip("xyz", (0.23 * a, 0.31 * b, 0.17 * d))]
    p1 = [snap(L[ax], v + cell) for ax, v in zip("xyz", p0)]
    exc = csx.AddExcitation("kick", exc_type=0, exc_val=[1, 1, 1])
    exc.AddBox(p0, p1)

    # probes: one-cell voltage lines along x, y, z at a different asymmetric spot
    q = [snap(L[ax], v) for ax, v in zip("xyz", (0.71 * a, 0.62 * b, 0.58 * d))]
    for i, ax in enumerate("xyz"):
        e = list(q); e[i] = snap(L[ax], q[i] + cell)
        csx.AddProbe(f"u{ax}", p_type=0).AddBox(q, e)

    path = sim_path(CASE)
    fdtd.Run(path, cleanup=True, verbose=0)

    spec = 0
    for ax in "xyz":
        t, u = np.loadtxt(Path(path) / f"u{ax}", comments="%", unpack=True)
        dt = t[1] - t[0]
        nfft = 1 << 22
        U = np.fft.rfft(u * np.hanning(len(u)), nfft)
        spec = spec + abs(U) / abs(U).max()
    freq = np.fft.rfftfreq(nfft, dt)
    band = (freq > f_lo * 0.97) & (freq < f_hi * 1.03)
    fb, sb = freq[band], 20 * np.log10(spec[band] + 1e-30)

    # local maxima above (max - 60 dB), parabolic refinement
    peaks = []
    thr = sb.max() - 60
    for i in range(1, len(sb) - 1):
        if sb[i] > sb[i - 1] and sb[i] >= sb[i + 1] and sb[i] > thr:
            y0, y1, y2 = sb[i - 1], sb[i], sb[i + 1]
            off = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
            peaks.append(fb[i] + off * (fb[1] - fb[0]))
    peaks = np.array(peaks)

    r = CaseResult(CASE, f"PEC cavity {a} x {b} x {d} mm, modes {f_lo/1e9:g}-{f_hi/1e9:g} GHz",
                   "Eigenmode solver, PEC boundaries, N modes")
    for f_an, names in analytic_modes():
        # nearest spectral peak; peaks must lie within 1% to count as found
        j = np.argmin(abs(peaks - f_an)) if len(peaks) else None
        f_num = peaks[j] if j is not None and abs(peaks[j] - f_an) / f_an < 0.01 else np.nan
        r.metrics.append(Metric(f"f({'/'.join(names)})", "GHz", f_num / 1e9, f_an / 1e9, "analytic", tol_pct=0.5))
    r.notes.append(f"{len(peaks)} spectral peaks found; NaN = mode not detected within 1% "
                   f"(T_sim = {t[-1]*1e9:.0f} ns, resolution ~{1/t[-1]/1e6:.0f} MHz)")
    add_hfss_metrics(r)
    r.runtime_s = time.time() - t0

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(fb / 1e9, sb, "k-", lw=1, label="openEMS probe spectrum")
    for f_an, names in analytic_modes():
        ax.axvline(f_an / 1e9, color="r", ls="--", lw=0.8)
        ax.text(f_an / 1e9, sb.max() + 2, "/".join(names), rotation=90, fontsize=7, va="bottom", ha="center")
    ax.set(xlabel="f (GHz)", ylabel="dB (norm.)", title="Cavity resonances: openEMS spectrum vs analytic (red)")
    ax.set_ylim(sb.max() - 80, sb.max() + 25); ax.grid()
    savefig(fig, CASE)
    return r


if __name__ == "__main__":
    r = main(); r.print(); r.save()
