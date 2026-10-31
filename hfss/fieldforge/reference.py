"""Reference results (e.g. exported from Ansys HFSS) and comparison metrics.

Supported files:

* Touchstone ``.sNp`` (v1/v2 options line; DB/MA/RI; any frequency unit)
* HFSS report CSV export, e.g. columns ``"Freq [GHz]","dB(S(1,1)) []"``
"""
from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass, field

import numpy as np

FREQ_UNITS = {"HZ": 1.0, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9, "THZ": 1e12}


@dataclass
class Trace:
    label: str
    freq: np.ndarray                 # Hz
    db: np.ndarray                   # magnitude in dB
    complex: np.ndarray | None = None


@dataclass
class Reference:
    path: str
    label: str
    traces: list[Trace] = field(default_factory=list)

    def reflection(self) -> Trace | None:
        """The S11-like trace (first diagonal entry)."""
        for t in self.traces:
            if re.search(r"S\(\s*([^,()]+)\s*,\s*\1\s*\)", t.label) or t.label.startswith("S11"):
                return t
        return self.traces[0] if self.traces else None


def load(path: str, label: str | None = None) -> Reference:
    ext = os.path.splitext(path)[1].lower()
    label = label or os.path.splitext(os.path.basename(path))[0]
    if re.fullmatch(r"\.s\d+p", ext) or ext == ".ts":
        return _load_touchstone(path, label)
    if ext in (".csv", ".txt"):
        return _load_csv(path, label)
    raise ValueError(f"Unsupported reference file type '{ext}' (use .sNp Touchstone or HFSS report .csv)")


def _load_touchstone(path: str, label: str) -> Reference:
    m = re.search(r"\.s(\d+)p$", path.lower())
    nports = int(m.group(1)) if m else None
    unit, fmt, z0 = 1e9, "MA", 50.0
    numbers: list[float] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.split("!", 1)[0].strip()
            if not line:
                continue
            if line.startswith("#"):
                toks = line[1:].upper().split()
                for i, t in enumerate(toks):
                    if t in FREQ_UNITS:
                        unit = FREQ_UNITS[t]
                    elif t in ("DB", "MA", "RI"):
                        fmt = t
                    elif t == "R" and i + 1 < len(toks):
                        z0 = float(toks[i + 1])
                    elif t in ("Y", "Z", "H", "G"):
                        raise ValueError(f"Touchstone parameter type '{t}' is not supported (need S)")
                continue
            if line.startswith("["):
                kw = line.lower()
                if kw.startswith("[number of ports]"):
                    nports = int(line.split("]")[1])
                continue
            numbers.extend(float(x) for x in line.split())
    if nports is None:
        raise ValueError("Cannot determine the number of ports (use a .sNp file name)")
    per = 1 + 2 * nports * nports
    data = np.array(numbers[: len(numbers) // per * per]).reshape(-1, per)
    f = data[:, 0] * unit
    pairs = data[:, 1:].reshape(len(f), nports * nports, 2)
    a, b = pairs[..., 0], pairs[..., 1]
    if fmt == "RI":
        s = a + 1j * b
    elif fmt == "MA":
        s = a * np.exp(1j * np.deg2rad(b))
    else:
        s = 10 ** (a / 20) * np.exp(1j * np.deg2rad(b))
    s = s.reshape(len(f), nports, nports)
    if nports == 2:
        s = s.transpose(0, 2, 1)       # 2-port files are ordered S11 S21 S12 S22
    ref = Reference(path, label)
    for i in range(nports):
        for j in range(nports):
            sij = s[:, i, j]
            ref.traces.append(Trace(f"S({i + 1},{j + 1})", f, 20 * np.log10(np.maximum(np.abs(sij), 1e-15)), sij))
    return ref


def _load_csv(path: str, label: str) -> Reference:
    with open(path, encoding="utf-8-sig", errors="replace") as fh:
        rows = [r for r in csv.reader(fh) if r]
    header, body = rows[0], rows[1:]
    data = np.array([[float(x) for x in r[: len(header)]] for r in body if _is_number(r[0])])
    fcol = 0
    m = re.search(r"\[(\w*Hz)\]", header[0], re.I)
    unit = FREQ_UNITS.get(m.group(1).upper(), 1e9) if m else 1e9
    f = data[:, fcol] * unit
    ref = Reference(path, label)
    for k, name in enumerate(header[1:], start=1):
        name = re.sub(r"\s*\[[^\]]*\]\s*$", "", name.strip())
        vals = data[:, k]
        if name.lower().startswith("db("):
            ref.traces.append(Trace(name[3:-1] if name.endswith(")") else name, f, vals))
        elif name.lower().startswith("mag("):
            ref.traces.append(Trace(name[4:-1], f, 20 * np.log10(np.maximum(vals, 1e-15))))
        else:
            ref.traces.append(Trace(name, f, vals))
    if not ref.traces:
        raise ValueError("No data columns found in the CSV file")
    return ref


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


# --------------------------------------------------------------------------- metrics

@dataclass
class Metrics:
    f_res: float            # Hz, frequency of the deepest dip
    depth_db: float
    bw_lo: float | None     # -10 dB band around the dip (Hz)
    bw_hi: float | None

    @property
    def bandwidth(self):
        return None if self.bw_lo is None else self.bw_hi - self.bw_lo


def metrics(freq: np.ndarray, db: np.ndarray, level: float = -10.0) -> Metrics:
    k = int(np.argmin(db))
    lo = hi = None
    if db[k] < level:
        i = k
        while i > 0 and db[i - 1] < level:
            i -= 1
        j = k
        while j < len(db) - 1 and db[j + 1] < level:
            j += 1
        lo = _cross(freq, db, i - 1, i, level) if i > 0 else freq[0]
        hi = _cross(freq, db, j, j + 1, level) if j < len(db) - 1 else freq[-1]
    return Metrics(float(freq[k]), float(db[k]), lo, hi)


def _cross(f, y, i, j, level):
    if y[j] == y[i]:
        return float(f[i])
    return float(f[i] + (level - y[i]) * (f[j] - f[i]) / (y[j] - y[i]))


def compare(f_a, db_a, f_b, db_b) -> dict:
    """Metrics for two |S11| curves plus their RMS difference on the common band."""
    ma, mb = metrics(f_a, db_a), metrics(f_b, db_b)
    lo, hi = max(f_a.min(), f_b.min()), min(f_a.max(), f_b.max())
    out = {"a": ma, "b": mb, "rms_db": None, "max_db": None, "band": (lo, hi)}
    if hi > lo:
        f = np.linspace(lo, hi, 801)
        d = np.interp(f, f_a, db_a) - np.interp(f, f_b, db_b)
        out["rms_db"] = float(np.sqrt(np.mean(d ** 2)))
        out["max_db"] = float(np.abs(d).max())
    out["df_res"] = ma.f_res - mb.f_res
    out["df_res_pct"] = 100 * (ma.f_res - mb.f_res) / mb.f_res
    return out
