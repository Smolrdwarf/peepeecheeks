"""Case 5 - adaptive mesh convergence (HFSS "Maximum Delta S") on the microstrip of case 3.

Starts from the openEMS-tutorial mesh, which is known to be too coarse (eps_eff +8.6 %),
and lets adaptive.py refine it until max|dS| < 0.02 between passes, as HFSS would with
its default convergence settings. Passes if the loop converges by itself AND the
converged result matches Hammerstad-Jensen, i.e. the delta-S criterion caught the error
without being told the answer.

Mesh scale s maps to: lambda/(50 s) cells, strip-edge cells res/(4 s), 5 s substrate cells.
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from common import CaseResult, Metric, add_hfss_metrics, savefig
from adaptive import adaptive_mesh
from cases import microstrip_line as msl

CASE = "microstrip_adaptive"
MAX_DELTA_S = 0.02


def mesh_kw(s):
    return dict(res_div=50 * s, edge_div=4 * s, nz_sub=int(round(5 * s)) + 1)


def simulate(s):
    out = msl.simulate(tag="_adaptive", **mesh_kw(s))
    out["S"] = np.stack([out["s11"], out["s21"]], axis=-1)   # symmetric 2-port: S22=S11, S12=S21
    return out


def main() -> CaseResult:
    t0 = time.time()
    res = adaptive_mesh(simulate, lambda s: msl.build(**mesh_kw(s))[2],
                        start_scale=1.0, refine_factor=1.25, max_delta=MAX_DELTA_S,
                        min_converged_passes=1, max_passes=6, max_cells=3_000_000)
    print(res.table())
    fin = res.final.output
    i = np.argmin(abs(fin["freq"] - msl.f_eval))
    z0_an, eeff_an = msl.hammerstad_jensen(msl.er, msl.w, msl.h)

    r = CaseResult(CASE, "Adaptive mesh loop (max|dS| < 0.02) on the case-3 microstrip, from the tutorial mesh",
                   "Adaptive passes with Maximum Delta S = 0.02, Minimum Converged Passes = 1")
    r.metrics += [
        Metric("final max|dS| between passes", "-", res.final.delta_s if res.final.delta_s is not None else np.inf,
               MAX_DELTA_S, "hfss-criterion", upper_bound=True),
        Metric("converged Z0 @ 2 GHz", "ohm", fin["z0"][i], z0_an, "analytic", tol_pct=2.0),
        Metric("converged eps_eff @ 2 GHz", "-", fin["eeff"][i], eeff_an, "analytic", tol_pct=2.0),
    ]
    r.notes.append(res.reason)
    for p in res.passes:
        j = np.argmin(abs(p.output["freq"] - msl.f_eval))
        r.notes.append(f"pass {p.index}: scale {p.scale:.3f}, {p.cells:,} cells, {p.runtime_s:.0f}s, "
                       f"max|dS| {'-' if p.delta_s is None else f'{p.delta_s:.4f}'}, "
                       f"Z0 {p.output['z0'][j]:.2f} ohm, eps_eff {p.output['eeff'][j]:.4f}")
    add_hfss_metrics(r)
    r.runtime_s = time.time() - t0

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    idx = [p.index for p in res.passes]
    ax[0].plot(idx[1:], [p.delta_s for p in res.passes[1:]], "o-"); ax[0].axhline(MAX_DELTA_S, color="r", ls="--", label="limit 0.02")
    ax[0].set_yscale("log"); ax[0].set(xlabel="pass", ylabel="max |dS|", title="Convergence"); ax[0].legend(); ax[0].grid(which="both")
    for k, (key, ref, lab) in enumerate((("eeff", eeff_an, "eps_eff"), ("z0", z0_an, "Z0 (ohm)")), start=1):
        ax[k].plot(idx, [p.output[key][np.argmin(abs(p.output["freq"] - msl.f_eval))] for p in res.passes], "o-", label="openEMS @ 2 GHz")
        ax[k].axhline(ref, color="k", ls="--", label="Hammerstad-Jensen"); ax[k].set(xlabel="pass", ylabel=lab, title=f"{lab} per pass")
        ax[k].legend(); ax[k].grid()
    savefig(fig, CASE)
    return r


if __name__ == "__main__":
    r = main(); r.print(); r.save()
