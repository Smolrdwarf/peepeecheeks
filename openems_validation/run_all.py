"""Run every validation case and write results/SUMMARY.md.

usage:  python run_all.py [case ...]      (default: all cases)
"""
import importlib
import json
import platform
import sys
import time
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS  # noqa: E402  (also sets the openEMS DLL path)

CASES = ["wr90_waveguide", "cavity_eigenmode", "microstrip_line", "pec_sphere_rcs", "microstrip_adaptive",
         "elmer_cavity_eigenmode", "elmer_wr90_waveguide"]


def fmt_tol(m):
    if m["upper_bound"]:
        return f"<= {m['reference']:.4g}"
    if m["tol_abs"] is not None:
        return f"±{m['tol_abs']:g} {m['unit']}"
    return f"±{m['tol_pct']:g} %"


def fmt_err(m):
    if m["tol_abs"] is not None or m["upper_bound"]:
        return f"{m['error_abs']:+.3g} {m['unit']}"
    return f"{m['error_pct']:+.3f} %"


def write_summary(names):
    import openEMS
    lines = [f"# openEMS validation summary",
             "",
             f"Generated {date.today()} by `run_all.py` — openEMS {getattr(openEMS, '__version__', '?')}, "
             f"Python {platform.python_version()}, {platform.system()} {platform.release()}",
             "",
             "| Case | Solver | Result | Runtime |", "|---|---|---|---|"]
    data = {}
    for n in names:
        f = RESULTS / f"{n}.json"
        if f.exists():
            data[n] = json.loads(f.read_text())
            lines.append(f"| `{n}` | {data[n].get('solver', 'openEMS')} | {'PASS' if data[n]['passed'] else '**FAIL**'} | {data[n]['runtime_s']:.0f} s |")
    for n, d in data.items():
        lines += ["", f"## `{n}` — {d['description']}", "",
                  f"HFSS equivalent: {d['hfss_equivalent']}", "",
                  f"![{n}]({n}.png)", "",
                  f"| Metric | {d.get('solver', 'openEMS')} | Reference | Source | Error | Tolerance | |",
                  "|---|---|---|---|---|---|---|"]
        for m in d["metrics"]:
            lines.append(f"| {m['name']} | {m['openems']:.5g} {m['unit']} | {m['reference']:.5g} | {m['reference_source']} "
                         f"| {fmt_err(m)} | {fmt_tol(m)} | {'ok' if m['passed'] else '**FAIL**'} |")
        lines += [""] + [f"- {x}" for x in d["notes"]]
    (RESULTS / "SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    names = sys.argv[1:] or CASES
    failed = []
    for n in names:
        t0 = time.time()
        print(f"--- running {n}", flush=True)
        try:
            r = importlib.import_module(f"cases.{n}").main()
            r.print(); r.save()
            if not r.passed:
                failed.append(n)
        except Exception as e:  # keep going so one broken case doesn't hide the others
            print(f"!!! {n} crashed after {time.time()-t0:.0f}s: {e!r}")
            failed.append(n)
    write_summary(CASES)
    print(f"\nwrote {RESULTS / 'SUMMARY.md'}; failed: {failed or 'none'}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
