#!/usr/bin/env python3
"""The report's tables, in Markdown, from out/eval.*.json and out/time.*.jsonl."""
import glob
import json
import os

W = os.path.dirname(os.path.abspath(__file__))
rows = {}
for f in glob.glob(os.path.join(W, "out", "eval.*.json")):
    d = json.load(open(f))
    rows[(d["label"], d["tol"])] = d

times = {}


def add_times(fname, keyfn):
    p = os.path.join(W, "out", fname)
    if not os.path.exists(p):
        return
    for line in open(p):
        line = line.strip()
        if not line.startswith("{"):
            continue
        d = json.loads(line)
        k = keyfn(d)
        if k:
            times[k] = d.get("ns_per_cubic", d.get("ns_per_cubic_incl_python_objects"))


add_times("time.rust.jsonl", lambda d: (d["method"], d["tol"]))
add_times("time.harfbuzz.jsonl", lambda d: ("harfbuzz", d["tol"]))
add_times("time.fonttools.jsonl", lambda d: ("fonttools_cu2qu", d["tol"]))
add_times("time.cubic2quad.jsonl", lambda d: ("cubic2quad_js", d["tol"]))
add_times("time.fontforge.jsonl", lambda d: ("fontforge", 1.0))
p = os.path.join(W, "out", "time.skia.jsonl")
if os.path.exists(p):
    for line in open(p):
        d = json.loads(line)
        t = d["tol"]
        for tol in (1.0, 0.5, 0.1):
            if abs(t - tol) < 1e-9:
                times[("skia", tol)] = d["ns_per_cubic"]
            if abs(t - 3.375 * tol) < 1e-6:
                times[("skia_cal", tol)] = d["ns_per_cubic"]

NAMES = [
    ("fonttools_cu2qu", "fontTools cu2qu (real, pure Python)"),
    ("kurbo_approx_spline", "kurbo approx_spline (real)"),
    ("harfbuzz", "HarfBuzz hb-gpu-cu2qu.hh (real header, shim)"),
    ("kurbo_to_quads", "kurbo CubicBez::to_quads (real)"),
    ("lyon", "lyon_geom 1.0.20 for_each_quadratic_bezier (real)"),
    ("lyon017", "lyon_geom 0.17.7 cubic_to_quadratics (real)"),
    ("lyon017_monotonic", "lyon_geom 0.17.7 cubic_to_monotonic_quadratics (real)"),
    ("skia", "Skia convertCubicToQuads, tolScale = tol (real)"),
    ("skia_cal", "Skia convertCubicToQuads, tolScale = 3.375 tol (real)"),
    ("cubic2quad_js", "cubic2quad JS (real)"),
    ("fontforge", "FontForge 20230101 (real; fixed tolerance 1 / 1.5)"),
    ("reimpl_truong_presub", "Truong pairs on m equal pieces, closed-form m (reimpl.)"),
    ("reimpl_colomitchi_half", "Colomitchi half-splitting (reimpl.)"),
    ("reimpl_colomitchi_adaptive", "Colomitchi simple adaptive (reimpl.)"),
    ("reimpl_colomitchi_symmetric", "Colomitchi symmetric adaptive (reimpl.)"),
    ("reimpl_alfke_penner", "Alfke / Penner tangent + midpoint metric (reimpl.)"),
    ("oracle_ramp", "oracle: cu2qu construction, least n by Hausdorff"),
    ("oracle_truong", "oracle: Truong pairs, least even n by Hausdorff"),
    ("oracle_qm", "oracle: mid-point construction, least n by Hausdorff"),
    ("oracle_tangent_implied", "oracle: tangent-intersection off-curves, implied on-curves"),
    ("oracle_tangent", "oracle: tangent-intersection, on-curves on the cubic"),
]

for tol in (1.0, 0.5, 0.1):
    ref = rows[("fonttools_cu2qu", tol)]["segs"]
    print(f"\n### tolerance {tol} font units ({rows[('fonttools_cu2qu', tol)]['cubics']} cubics)\n")
    print("| method | quads per cubic | vs cu2qu | TT points per cubic | Hausdorff/tol mean | p99 | max | cubics over tol | end tangent off > 1 deg | joint kink > 1 deg | ns per cubic |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for key, label in NAMES:
        d = rows.get((key, tol))
        if not d:
            continue
        t = times.get((key, tol))
        ts = f"{t:,.0f}" if t else ""
        n_tan = d["cubics"] - d["failed"]
        print(
            f"| {label} | {d['mean_n']:.3f} | {d['segs'] / ref:.3f} | {d['tt_mean']:.3f} | {d['haus']['mean']:.3f} | {d['haus']['p99']:.3f} | "
            f"{d['haus']['max']:.3f} | {d['viol']} | {d['tan_gt_1']} ({100 * d['tan_gt_1'] / n_tan:.1f}%) | {d['kink_gt_1']} | {ts} |"
        )

print("\n### distribution of quadratics per cubic (share of cubics, %)\n")
print("| method | tol | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | >8 |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for key, label in NAMES:
    for tol in (1.0, 0.5, 0.1):
        d = rows.get((key, tol))
        if not d or key.startswith("reimpl_colo") or key in ("lyon017_monotonic", "reimpl_alfke_penner"):
            continue
        tot = sum(d["hist"].values())
        cells = []
        for k in range(1, 9):
            cells.append(f"{100 * d['hist'].get(str(k), 0) / tot:.1f}")
        more = sum(v for k, v in d["hist"].items() if int(k) > 8)
        cells.append(f"{100 * more / tot:.1f}")
        print(f"| {label} | {tol} | " + " | ".join(cells) + " |")
