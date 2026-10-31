"""Shared helpers for the openEMS-vs-HFSS validation suite."""
from __future__ import annotations

import json
import os
import sys
import tempfile
from dataclasses import dataclass, field, asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
HFSS_REF = ROOT / "references" / "hfss"

# openEMS Windows build: the Python modules need to find CSXCAD.dll / openEMS.dll.
_DEFAULT_INSTALL = Path.home() / "tools" / "openEMS"
for var in ("CSXCAD_INSTALL_PATH", "OPENEMS_INSTALL_PATH"):
    os.environ.setdefault(var, str(_DEFAULT_INSTALL))


def sim_path(case: str) -> str:
    """Scratch directory for raw solver output (kept out of the repo)."""
    p = Path(tempfile.gettempdir()) / "openems_validation" / case
    p.mkdir(parents=True, exist_ok=True)
    return str(p)


@dataclass
class Metric:
    name: str
    unit: str
    openems: float
    reference: float
    reference_source: str  # "analytic" | "hfss" | "literature"
    tol_pct: float | None = None   # relative tolerance in percent
    tol_abs: float | None = None   # absolute tolerance in `unit`
    upper_bound: bool = False      # if True, reference is a limit: pass when openems <= reference

    def __post_init__(self):
        self.openems = float(self.openems)
        self.reference = float(self.reference)

    @property
    def error_abs(self) -> float:
        return self.openems - self.reference

    @property
    def error_pct(self) -> float:
        return 100.0 * self.error_abs / self.reference if self.reference else float("nan")

    @property
    def passed(self) -> bool:
        if self.upper_bound:
            return bool(self.openems <= self.reference)
        if self.tol_abs is not None:
            return bool(abs(self.error_abs) <= self.tol_abs)
        if self.tol_pct is not None:
            return bool(abs(self.error_pct) <= self.tol_pct)
        return True


@dataclass
class CaseResult:
    case: str
    description: str
    hfss_equivalent: str
    metrics: list[Metric] = field(default_factory=list)
    runtime_s: float = 0.0
    notes: list[str] = field(default_factory=list)
    solver: str = "openEMS"        # which open-source solver produced Metric.openems

    @property
    def passed(self) -> bool:
        return all(m.passed for m in self.metrics)

    def save(self) -> Path:
        RESULTS.mkdir(exist_ok=True)
        out = RESULTS / f"{self.case}.json"
        d = asdict(self)
        for md, m in zip(d["metrics"], self.metrics):
            md.update(error_abs=m.error_abs, error_pct=m.error_pct, passed=m.passed)
        d["passed"] = self.passed
        out.write_text(json.dumps(d, indent=2))
        return out

    def print(self) -> None:
        print(f"\n=== {self.case} [{self.solver}]: {'PASS' if self.passed else 'FAIL'}  ({self.runtime_s:.1f} s)")
        for m in self.metrics:
            tol = ("<= ref" if m.upper_bound else
                   f"+/-{m.tol_abs} {m.unit}" if m.tol_abs is not None else f"+/-{m.tol_pct}%")
            err = (f"{m.error_abs:+.3g} {m.unit}" if (m.tol_abs is not None or m.upper_bound)
                   else f"{m.error_pct:+.3f}%")
            print(f"  {m.name:<34} {self.solver}={m.openems:12.5g}  ref({m.reference_source})={m.reference:12.5g} "
                  f"{m.unit:<6} err={err}  tol {tol}  {'ok' if m.passed else 'FAIL'}")
        for n in self.notes:
            print("  note:", n)


def load_hfss_csv(case: str) -> dict[str, float] | None:
    """Optional HFSS reference values: references/hfss/<case>.csv with `metric,value` rows.

    Metric names must match the `name` field used by the case script.
    """
    f = HFSS_REF / f"{case}.csv"
    if not f.exists():
        return None
    vals = {}
    for line in f.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.lower().startswith("metric,"):
            continue
        k, v = line.rsplit(",", 1)
        vals[k.strip()] = float(v)
    return vals


def add_hfss_metrics(res: CaseResult) -> None:
    """Duplicate each metric against HFSS reference values if they have been supplied."""
    ref = load_hfss_csv(res.case)
    if ref is None:
        res.notes.append(f"no HFSS reference supplied (references/hfss/{res.case}.csv); compared to analytic only")
        return
    extra = []
    for m in res.metrics:
        if m.name in ref:
            extra.append(Metric(m.name, m.unit, m.openems, ref[m.name], "hfss", m.tol_pct, m.tol_abs, m.upper_bound))
    res.metrics.extend(extra)


def savefig(fig, case: str, suffix: str = "") -> None:
    RESULTS.mkdir(exist_ok=True)
    fig.savefig(RESULTS / f"{case}{suffix}.png", dpi=110, bbox_inches="tight")
