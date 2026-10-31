"""Unit handling for AEDT expressions.

AEDT expressions mix numbers with unit suffixes and variables, e.g. ``'2*a + 3mm'``,
``'90deg'``, ``'2.87GHz'``. FieldForge's evaluator is unit-less, so unit-bearing
literals are rewritten into canonical units before evaluation:

* lengths      -> the design's model units (e.g. mm)
* angles       -> degrees
* frequencies  -> GHz
* everything else (ohm, S/m, ...) -> SI value
"""
from __future__ import annotations

import math
import re

LENGTH = {"fm": 1e-15, "pm": 1e-12, "nm": 1e-9, "um": 1e-6, "mm": 1e-3, "cm": 1e-2, "dm": 1e-1,
          "meter": 1.0, "m": 1.0, "km": 1e3, "uin": 25.4e-9, "mil": 25.4e-6, "in": 25.4e-3,
          "inch": 25.4e-3, "ft": 0.3048, "yd": 0.9144, "mile": 1609.344}
ANGLE = {"deg": 1.0, "degree": 1.0, "rad": 180.0 / math.pi, "mrad": 0.18 / math.pi}
FREQ = {"Hz": 1e-9, "kHz": 1e-6, "MHz": 1e-3, "GHz": 1.0, "THz": 1e3}
OTHER = {"ohm": 1.0, "Ohm": 1.0, "kOhm": 1e3, "kohm": 1e3, "mOhm": 1e-3, "MOhm": 1e6,
         "siemens/m": 1.0, "S_per_meter": 1.0, "cel": 1.0, "kel": 1.0, "W": 1.0, "mW": 1e-3,
         "V": 1.0, "mV": 1e-3, "A": 1.0, "mA": 1e-3, "s": 1.0, "ns": 1e-9, "ps": 1e-12,
         "H": 1.0, "nH": 1e-9, "pH": 1e-12, "F": 1.0, "pF": 1e-12, "nF": 1e-9, "uF": 1e-6}

_ALL_UNITS = sorted({*LENGTH, *ANGLE, *FREQ, *OTHER}, key=len, reverse=True)
_NUM_UNIT = re.compile(
    r"(?<![\w.$])((?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*(" + "|".join(map(re.escape, _ALL_UNITS)) + r")(?!\w)")
_PAREN_UNIT = re.compile(r"\)\s*(" + "|".join(map(re.escape, _ALL_UNITS)) + r")(?!\w)")


def unit_scale(unit: str, model_units: str) -> float:
    if unit in LENGTH:
        return LENGTH[unit] / LENGTH[model_units]
    if unit in ANGLE:
        return ANGLE[unit]
    if unit in FREQ:
        return FREQ[unit]
    return OTHER.get(unit, 1.0)


def var_name(name: str) -> str:
    """AEDT project variables start with '$', which is not a valid identifier."""
    return "prj_" + name[1:] if name.startswith("$") else name


def convert_expr(expr: str, model_units: str) -> str:
    """Rewrite an AEDT expression into a unit-less one in canonical units."""
    s = str(expr).strip()
    s = s.replace("^", "**")
    s = re.sub(r"\$(\w+)", r"prj_\1", s)
    s = re.sub(r"\bif\s*\(", "_if(", s)  # 'if' is a Python keyword

    def repl(m):
        val = float(m.group(1)) * unit_scale(m.group(2), model_units)
        return f"({val!r})"

    s = _NUM_UNIT.sub(repl, s)
    # a unit applied to a parenthesised group: '(0.25*2) mm'
    return _PAREN_UNIT.sub(lambda m: f")*({unit_scale(m.group(1), model_units)!r})", s)


def split_value_unit(text: str) -> tuple[float, str] | None:
    """'2.87GHz' -> (2.87, 'GHz'); returns None for non-literal expressions."""
    m = re.fullmatch(r"\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*([A-Za-z_/]*)\s*", str(text))
    if not m:
        return None
    return float(m.group(1)), m.group(2)
