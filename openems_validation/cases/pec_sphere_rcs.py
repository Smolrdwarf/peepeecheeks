"""Case 4 - monostatic radar cross section of a PEC sphere (radius 200 mm), 50 MHz-1 GHz.

HFSS equivalent: Driven Modal with incident plane wave, radiation boundary or PML box,
far-field monostatic RCS sweep. Reference: exact Mie series.

openEMS: total-field/scattered-field plane wave box, NF2FF box, backscatter direction.
Note the sphere is staircased on the Cartesian FDTD grid (HFSS uses curved tetrahedra).

Mesh study (2026-10-31): uniform lambda/25 (12 mm) gives 0.83 dB RMS error; a 6 mm
region around the sphere gives 0.40 dB. Remaining error is largest next to the deep
nulls (ka ~ 1.75, 2.95), where a small shift in ka produces a large dB difference.
A 4 mm run was attempted but stopped (host low on memory), so convergence below 6 mm is unverified.
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.special import spherical_jn, spherical_yn
from common import CaseResult, Metric, sim_path, add_hfss_metrics, savefig

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0, Z0
from openEMS.ports import UI_data

CASE = "pec_sphere_rcs"
unit = 1e-3
R = 200.0                     # sphere radius, mm
f_start, f_stop = 50e6, 1000e6
ka_check = [0.5, 1.0, 2.0, 3.0]


def mie_monostatic(ka, nmax=None):
    """Normalized monostatic RCS sigma/(pi a^2) of a PEC sphere."""
    ka = np.atleast_1d(ka).astype(float)
    out = np.empty_like(ka)
    for i, x in enumerate(ka):
        N = nmax or int(x + 4 * x ** (1 / 3) + 10)
        n = np.arange(1, N + 1)
        jn, yn = spherical_jn(n, x), spherical_yn(n, x)
        djn, dyn = spherical_jn(n, x, True), spherical_yn(n, x, True)
        hn, dhn = jn - 1j * yn, djn - 1j * dyn
        a = jn / hn
        b = (jn + x * djn) / (hn + x * dhn)
        s = np.sum((-1) ** n * (2 * n + 1) * (a - b))
        out[i] = abs(s) ** 2 / x ** 2           # sigma/(pi a^2) = |S|^2 / (ka)^2
    return out


def main(cells_per_lambda=25, fine_cell=6.0) -> CaseResult:
    t0 = time.time()
    sim_box, pw_box = 1200.0, 750.0
    fdtd = openEMS(EndCriteria=1e-5)
    fdtd.SetGaussExcite(0.5 * (f_start + f_stop), 0.5 * (f_stop - f_start))
    fdtd.SetBoundaryCond(['PML_8'] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(unit)
    if fine_cell:   # fine uniform region around the sphere, graded to coarse outside
        nf = int(np.ceil(1.15 * R / fine_cell))
        mesh.SetLines('x', np.r_[np.linspace(-nf, nf, 2 * nf + 1) * fine_cell, -sim_box / 2, sim_box / 2])
    else:
        mesh.SetLines('x', [-sim_box / 2, 0, sim_box / 2])
    mesh.SmoothMeshLines('x', C0 / f_stop / unit / cells_per_lambda, ratio=1.3)
    mesh.SetLines('y', mesh.GetLines('x')); mesh.SetLines('z', mesh.GetLines('x'))

    csx.AddMetal('sphere').AddSphere(priority=10, center=[0, 0, 0], radius=R)
    pw = csx.AddExcitation('plane_wave', exc_type=10, exc_val=[0, 0, 1])   # E_z, travelling +x
    pw.SetPropagationDir([1, 0, 0])
    pw.SetFrequency(0.5 * (f_start + f_stop))
    pw.AddBox([-pw_box / 2] * 3, [pw_box / 2] * 3)
    nf2ff = fdtd.CreateNF2FFBox()

    path = sim_path(CASE)
    fdtd.Run(path, cleanup=True, verbose=0)

    freq = np.linspace(f_start, f_stop, 191)
    ef = UI_data('et', path, freq)
    Pin = 0.5 / Z0 * abs(np.array(ef.ui_f_val[0])) ** 2
    res = nf2ff.CalcNF2FF(path, freq, 90, 180, outfile='back_nf2ff.h5', read_cached=False)
    rcs = np.array([4 * np.pi / Pin[i] * res.P_rad[i][0][0] for i in range(len(freq))])
    ka = 2 * np.pi * freq / C0 * R * unit
    norm_num = rcs / (np.pi * (R * unit) ** 2)
    norm_mie = mie_monostatic(ka)
    db = lambda v: 10 * np.log10(v)

    r = CaseResult(CASE, f"PEC sphere a={R:g} mm, monostatic RCS, ka {ka[0]:.2f}-{ka[-1]:.2f}",
                   "Driven Modal, incident plane wave, monostatic RCS far-field sweep")
    for kc in ka_check:
        i = np.argmin(abs(ka - kc))
        r.metrics.append(Metric(f"RCS/(pi a^2) @ ka={kc:g}", "dB", db(norm_num[i]), db(norm_mie[i]), "analytic", tol_abs=0.5))
    band = ka >= 0.3
    rms = np.sqrt(np.mean((db(norm_num[band]) - db(norm_mie[band])) ** 2))
    r.metrics.append(Metric("RMS error, ka 0.3-4.2", "dB", rms, 0.5, "analytic", upper_bound=True))
    r.notes.append(f"FDTD cell = lambda/{cells_per_lambda} at {f_stop/1e9:g} GHz"
                   + (f", {fine_cell:g} mm around sphere" if fine_cell else "") + "; staircased sphere")
    add_hfss_metrics(r)
    r.runtime_s = time.time() - t0

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    kk = np.linspace(ka[0], ka[-1], 600)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.semilogy(kk, mie_monostatic(kk), "k-", label="Mie series (exact)")
    ax.semilogy(ka, norm_num, "r--", label="openEMS")
    ax.set(xlabel="ka", ylabel="sigma / (pi a^2)", title="PEC sphere monostatic RCS"); ax.legend(); ax.grid(which="both")
    savefig(fig, CASE)
    return r


if __name__ == "__main__":
    r = main(); r.print(); r.save()
