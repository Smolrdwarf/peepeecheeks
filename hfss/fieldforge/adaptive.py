"""HFSS-style adaptive mesh convergence for openEMS.

HFSS refines its tetrahedral mesh pass by pass and stops once the S-parameters change by
less than "Maximum Delta S" (default 0.02) between consecutive passes, for a given number
of "Minimum Converged Passes" (default 1). openEMS (FDTD) has no adaptive meshing, so we
emulate the loop from the outside: rebuild the model with a denser mesh at each pass and
compare S-parameters with the previous pass.

Differences from HFSS, deliberately:
  * Delta S is taken as the max over the whole frequency band (HFSS: one adaptive
    frequency). Stricter, and natural for FDTD, which is broadband anyway.
  * Refinement is global (one density scale), not error-driven local refinement.
  * Every pass is a full rerun; FDTD cannot reuse the previous solution.

The model is described by two callables so any geometry can use it:
    estimate_cells(scale) -> int                      # cheap, no solve
    simulate(scale)       -> dict with key "S": complex array (n_freq, ...) plus anything else
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable

import numpy as np


@dataclass
class Pass:
    index: int
    scale: float
    cells: int
    runtime_s: float
    delta_s: float | None          # None on the first pass
    output: dict = field(repr=False, default_factory=dict)


@dataclass
class AdaptiveResult:
    passes: list[Pass]
    converged: bool
    reason: str

    @property
    def final(self) -> Pass:
        return self.passes[-1]

    def table(self) -> str:
        rows = ["pass  scale   cells       time    max|dS|"]
        for p in self.passes:
            ds = "-" if p.delta_s is None else f"{p.delta_s:.4f}"
            rows.append(f"{p.index:>4}  {p.scale:5.2f}  {p.cells:>9,}  {p.runtime_s:6.1f}s  {ds}")
        rows.append(f"{'converged' if self.converged else 'NOT converged'}: {self.reason}")
        return "\n".join(rows)


def max_delta_s(s_new: np.ndarray, s_old: np.ndarray) -> float:
    """HFSS definition: max over all S-matrix entries (and here, all frequencies) of |S_n - S_n-1|."""
    return float(np.max(np.abs(np.asarray(s_new) - np.asarray(s_old))))


def adaptive_mesh(simulate: Callable[[float], dict],
                  estimate_cells: Callable[[float], int] | None = None,
                  *,
                  start_scale: float = 1.0,
                  refine_factor: float = 1.25,
                  max_delta: float = 0.02,
                  min_converged_passes: int = 1,
                  max_passes: int = 6,
                  max_cells: int | None = 5_000_000,
                  log: Callable[[str], None] = print) -> AdaptiveResult:
    """Run passes at scale = start_scale * refine_factor**k until max|dS| < max_delta
    for `min_converged_passes` consecutive passes, or a pass/cell limit is hit."""
    if refine_factor <= 1:
        raise ValueError("refine_factor must be > 1")
    passes: list[Pass] = []
    streak = 0
    scale = start_scale
    for k in range(1, max_passes + 1):
        if estimate_cells is not None and max_cells is not None:
            n = estimate_cells(scale)
            if n > max_cells:
                reason = f"next pass would need {n:,} cells > max_cells {max_cells:,}"
                log(f"[adaptive] stop: {reason}")
                return AdaptiveResult(passes, False, reason)
        t0 = time.time()
        out = simulate(scale)
        if "S" not in out:
            raise KeyError("simulate() must return a dict containing 'S'")
        ds = max_delta_s(out["S"], passes[-1].output["S"]) if passes else None
        cells = int(out.get("n_cells", estimate_cells(scale) if estimate_cells else 0))
        passes.append(Pass(k, scale, cells, time.time() - t0, ds, out))
        log(f"[adaptive] pass {k}: scale {scale:.3f}, {cells:,} cells, "
            f"{passes[-1].runtime_s:.0f}s, max|dS| = {'-' if ds is None else f'{ds:.4f}'}")
        if ds is not None and ds < max_delta:
            streak += 1
            if streak >= min_converged_passes:
                return AdaptiveResult(passes, True,
                                      f"max|dS| {ds:.4f} < {max_delta} for {streak} consecutive pass(es)")
        else:
            streak = 0
        scale *= refine_factor
    return AdaptiveResult(passes, False, f"max_passes ({max_passes}) reached")
