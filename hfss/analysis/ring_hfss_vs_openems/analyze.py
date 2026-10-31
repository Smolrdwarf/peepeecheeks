"""Tabulate every openEMS variant against the HFSS solved points and draw the attribution figure.

usage: python analyze.py   (reads hfss_extracted.json and runs/*/s11.npz, writes results.md + figure)
"""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
h = json.load(open(os.path.join(HERE, "hfss_extracted.json")))
pts = h["DV1666_SOL415_M0_V2505.sd"]["points"]
fh = np.array([p["f"] for p in pts])
sh = np.array([complex(p["re"], p["im"]) for p in pts])
interior = slice(1, -1)          # the 8 points away from the band edges


def at(f, s, fq):
    return np.interp(fq, f, s.real) + 1j * np.interp(fq, f, s.imag)


def unwrap_span(f, s):
    ph = np.unwrap(np.angle(s))
    return np.degrees(ph[-1] - ph[0])


def ceq(f0, s):
    z = 50 * (1 + s) / (1 - s)
    return -1 / (2 * np.pi * f0 * z.imag) * 1e12 if z.imag < 0 else float("nan")


rows = []
data = {}
for d in sorted(glob.glob(os.path.join(HERE, "runs", "*", "s11.npz"))):
    name = os.path.basename(os.path.dirname(d))
    r = np.load(d)
    meta = json.load(open(os.path.join(os.path.dirname(d), "meta.json")))
    f, s = r["f"], r["S11"]
    sq = at(f, s, fh)
    dist_int = np.abs(sq[interior] - sh[interior]).mean()
    dist_edge = np.abs(sq[[0, -1]] - sh[[0, -1]])
    data[name] = (f, s)
    rows.append({
        "name": name, "knobs": meta["knobs"], "cells": meta["cells"], "min_cell": meta["min_cell_mm"],
        "runtime": meta["runtime_s"], "dist_int": dist_int, "dist_lo": dist_edge[0], "dist_hi": dist_edge[1],
        "rot": unwrap_span(fh[interior], sq[interior]),
        "ceq": ceq(fh[3], sq[3]),
        "db_mid": 20 * np.log10(np.abs(sq[interior])).mean(),
    })

hfss_row = {"rot": unwrap_span(fh[interior], sh[interior]), "ceq": ceq(fh[3], sh[3]),
            "db_mid": 20 * np.log10(np.abs(sh[interior])).mean()}

lines = ["| variant | knob | cells | min cell (mm) | mean |dS11| interior | |dS11| @1.435 / @4.305 GHz | "
         "phase rotation 1.61-3.95 GHz (deg) | C_eq @2.15 GHz (pF) | mean |S11| interior (dB) |",
         "|---|---|---|---|---|---|---|---|---|",
         f"| **HFSS** | - | 16,155 tets | 0.115 (edge) | 0 | 0 / 0 | {hfss_row['rot']:.0f} | {hfss_row['ceq']:.3f} | "
         f"{hfss_row['db_mid']:.2f} |"]
for r in sorted(rows, key=lambda r: r["dist_int"]):
    knob = ", ".join(f"{k}={v}" for k, v in r["knobs"].items()) or "(none)"
    lines.append(f"| {r['name']} | {knob} | {r['cells']:,} | {r['min_cell']:.3f} | {r['dist_int']:.3f} | "
                 f"{r['dist_lo']:.2f} / {r['dist_hi']:.2f} | {r['rot']:.0f} | {r['ceq']:.3f} | {r['db_mid']:.2f} |")
table = "\n".join(lines)
open(os.path.join(HERE, "variants_table.md"), "w").write(table + "\n")
print(table)

# ------------------------------------------------------------------ figure
# survivor: best variant with the port gap open; runner-up: best variant with the gap closed
best = min((r for r in rows if "gap" in r["knobs"]), key=lambda r: r["dist_int"])["name"]
closed = [r for r in rows if "gap" not in r["knobs"] and r["name"] != "baseline"]
runner = min(closed, key=lambda r: r["dist_int"])["name"] if closed else None
fig = plt.figure(figsize=(12, 5.2), layout="constrained")
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
t = np.linspace(0, 2 * np.pi, 400)
ax2.plot(np.cos(t), np.sin(t), color="#888", lw=1)
for rr in (0.2, 1, 5):
    ax2.plot(rr / (1 + rr) + np.cos(t) / (1 + rr), np.sin(t) / (1 + rr), color="#ddd", lw=0.6)
for x in (0.5, 1, 2):
    for sx in (x, -x):
        z = 1j * sx + np.logspace(-3, 3, 400)
        g = (z - 1) / (z + 1)
        ax2.plot(g.real, g.imag, color="#ddd", lw=0.6)
style = {"baseline": ("#d62728", "-", "FieldForge baseline (gap snapped shut)"),
         best: ("#1f77b4", "-", f"FieldForge, {best} (HFSS port gap kept open)"),
         runner: ("#2ca02c", ":", f"FieldForge, {runner} (gap closed, finest mesh: fails)")}
for name, (col, ls, lab) in style.items():
    if name is None or name not in data:
        continue
    f, s = data[name]
    ax1.plot(f / 1e9, 20 * np.log10(np.abs(s)), color=col, ls=ls, lw=1.6, label=lab)
    ax2.plot(s.real, s.imag, color=col, ls=ls, lw=1.6)
ax1.plot(fh / 1e9, 20 * np.log10(np.abs(sh)), "ko", ms=5, label="HFSS solved points")
ax2.plot(sh.real, sh.imag, "ko", ms=5)
for i in (0, len(fh) - 1):
    ax2.annotate(f"{fh[i] / 1e9:.2f} GHz", (sh[i].real, sh[i].imag), fontsize=7, xytext=(4, 4),
                 textcoords="offset points")
ax1.set_xlabel("Frequency (GHz)")
ax1.set_ylabel("|S11| (dB)")
ax1.grid(True, color="#e0e0e0")
ax1.legend(fontsize=8, loc="lower left")
ax1.set_title("|S11|: HFSS vs openEMS variants", fontsize=10, loc="left")
ax2.set_aspect("equal")
ax2.axis("off")
ax2.set_title("Smith chart (same traces)", fontsize=10, loc="left")
fig.savefig(os.path.join(HERE, "attribution.png"), dpi=130)
print("best:", best, "runner-up:", runner)
