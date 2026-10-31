"""Extract the S11 points HFSS actually solved from a raw .aedtresults folder.

HFSS writes the interpolating-sweep solution as a text tree (DV*_SOL*_M0_*.sd). Each
directly solved frequency appears as a 'freq' block with its S-matrix. Adaptive-pass
solutions (SOL399 = last adaptive pass, SOL400 = adaptive solution) and the
convergence history (*.EMConvData) are extracted too.

usage: python extract_hfss.py <HFSSDesign1.results dir> <out.json>
"""
import glob
import json
import os
import re
import sys

C = re.compile(r"c\(([-\d.eE+]+),\s*([-\d.eE+]+)\)")


def solved_points(path):
    text = open(path, encoding="latin-1").read()
    pts = []
    # top-level 'freq' blocks and the 'DataPoint' blocks inside 'soln'
    for blk in re.split(r"\$begin '(?:freq|DataPoint)'", text)[1:]:
        f = re.search(r"Frequency=([\d.eE+]+)", blk)
        s = re.search(r"\$begin 'SMatrix'(.*?)\$end 'SMatrix'", blk, re.S)
        if not (f and s):
            continue
        re_, im = map(float, C.search(s.group(1)).groups())
        zw = re.search(r"\$begin 'Zwave'(.*?)\$end", blk, re.S)
        pts.append({"f": float(f.group(1)), "re": re_, "im": im,
                    "zwave": list(map(float, C.search(zw.group(1)).groups())) if zw else None})
    meta = {}
    for k in ("SolutionType", "BasisSize", "NumSamples", "ExtrapolationIsConverged", "ErrorTolerance",
              "FreqSequence", "AdaptivePass"):
        m = re.search(rf"{k}[=(]([^\n]*)", text)
        if m:
            meta[k] = m.group(1).strip().rstrip(")")
    uniq = {round(p["f"]): p for p in pts}
    return [uniq[k] for k in sorted(uniq)], meta


def convergence(path):
    text = open(path, encoding="latin-1").read()
    out = {"IsConverged": re.search(r"IsConverged=(\w+)", text).group(1),
           "ConvTarget": float(re.search(r"ConvTarget='([\d.]+)'", text).group(1)), "passes": []}
    for blk in re.findall(r"\$begin 'QuantityValues'(.*?)\$end 'QuantityValues'", text, re.S):
        if "Intrinsics=''" not in blk:
            continue
        n = int(re.search(r"PassNo=(\d+)", blk).group(1))
        tets = float(re.search(r"SolvedElements\(([\d.eE+]+)", blk).group(1))
        d = re.search(r"MaxMagDeltaS\(([\d.eE+-]+)", blk)
        out["passes"].append({"pass": n, "tets": tets, "maxMagDeltaS": float(d.group(1)) if d else None})
    return out


def main(results_dir, out):
    data = {}
    for sd in sorted(glob.glob(os.path.join(results_dir, "DV*_SOL*_M*_V*.sd"))):
        pts, meta = solved_points(sd)
        data[os.path.basename(sd)] = {"meta": meta, "points": pts}
    conv = glob.glob(os.path.join(results_dir, "*.EMConvData"))
    data["convergence"] = convergence(conv[0]) if conv else None
    json.dump(data, open(out, "w"), indent=1)
    for k, v in data.items():
        if k == "convergence":
            print("convergence:", v)
            continue
        print(k, v["meta"])
        for p in v["points"]:
            import math
            mag = math.hypot(p["re"], p["im"])
            print(f"   {p['f'] / 1e9:8.4f} GHz  S11 = {p['re']:+.4f}{p['im']:+.4f}j  |S11| = {20 * math.log10(mag):6.2f} dB")


if __name__ == "__main__":
    main(*sys.argv[1:3])
