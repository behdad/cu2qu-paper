#!/usr/bin/env python3
"""Questions answered from the evaluator's per-cubic output (B/percubic/*.txt)."""
import json
import os
import sys

import numpy as np

S = "/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad"
B = "/home/behdad/hb/cu2qu-work/research-survey"

cubics = np.array(json.load(open(os.path.join(S, "ft_cubics.json"))), dtype=float)
N = len(cubics)
zero_start = (cubics[:, 0] == cubics[:, 1]).all(1)
zero_end = (cubics[:, 2] == cubics[:, 3]).all(1)
zero_handle = zero_start | zero_end
d3 = np.hypot(*(cubics[:, 3] - 3 * cubics[:, 2] + 3 * cubics[:, 1] - cubics[:, 0]).T)
# collinear control polygon (a straight "curve")
v1 = cubics[:, 1] - cubics[:, 0]
v2 = cubics[:, 2] - cubics[:, 0]
v3 = cubics[:, 3] - cubics[:, 0]
cross = lambda a, b: a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
straight = (cross(v1, v3) == 0) & (cross(v2, v3) == 0)


def load(name):
    a = np.loadtxt(os.path.join(B, "percubic", name + ".txt"), skiprows=1)
    cols = "idx nseg nquad tt haus_c2s haus_s2c eq tan_start tan_end kink bbox d3".split()
    return {c: a[:, i] for i, c in enumerate(cols)}


def section(t):
    print("\n" + t)
    print("-" * len(t))


section(f"corpus: {N} cubics")
print(f"  zero-length handle at an end (p1 == p0 or p2 == p3): {zero_handle.sum()} ({100 * zero_handle.mean():.2f}%)")
print(f"  straight (all four points collinear):                 {straight.sum()} ({100 * straight.mean():.2f}%)")
print(f"  |d3| percentiles 50/90/99/max: {np.percentile(d3, 50):.1f} {np.percentile(d3, 90):.1f} {np.percentile(d3, 99):.1f} {d3.max():.1f}")

section("end tangents: which cubics lose them (angle between the cubic's and the spline's end tangent > 1 degree)")
for name in (
    "fonttools_cu2qu",
    "harfbuzz",
    "kurbo_to_quads",
    "cubic2quad_js",
    "skia",
    "skia_cal",
    "fontforge",
    "reimpl_truong_presub",
    "oracle_tangent_implied",
    "oracle_tangent",
):
    for tol in ("1.0", "0.5", "0.1"):
        p = os.path.join(B, "percubic", f"{name}.{tol}.txt")
        if not os.path.exists(p):
            continue
        d = load(f"{name}.{tol}")
        tan = np.maximum(d["tan_start"], d["tan_end"])
        bad = tan > 1.0
        n1 = d["nseg"] == 1
        print(
            f"  {name:22s} tol {tol}: >1deg {bad.sum():6d}  of which single-quad {int((bad & n1).sum()):6d}, "
            f"zero-handle cubic {int((bad & zero_handle).sum()):6d}, neither {int((bad & ~n1 & ~zero_handle).sum()):6d};  "
            f">10deg {int((tan > 10).sum()):6d}  >90deg {int((tan > 90).sum()):5d};  "
            f"multi-quad & nonzero handles: max {tan[~n1 & ~zero_handle].max() if (~n1 & ~zero_handle).any() else 0:.2e} deg"
        )

section("FontForge against its own thresholds (Hausdorff distance, font units)")
d = load("fontforge.1.0")
h = np.maximum(d["haus_c2s"], d["haus_s2c"])
for thr in (1.0, 1.5, 2.0, 3.0, 4.5):
    print(f"  > {thr}: {int((h > thr).sum()):6d} cubics ({100 * (h > thr).mean():.2f}%)")
print(f"  max {h.max():.3f}; mean {h.mean():.3f}; median {np.median(h):.3f}")

section("true (Hausdorff) distance over equal-parameter distance, uniform-split methods")
for name in ("fonttools_cu2qu", "kurbo_to_quads", "reimpl_truong_presub"):
    for tol in ("1.0", "0.5", "0.1"):
        d = load(f"{name}.{tol}")
        h = np.maximum(d["haus_c2s"], d["haus_s2c"])
        eq = d["eq"]
        ok = eq > 1e-9
        r = h[ok] / eq[ok]
        multi = ok & (d["nseg"] > 1)
        rm = h[multi] / eq[multi]
        print(
            f"  {name:22s} tol {tol}: H/eq mean {r.mean():.3f} median {np.median(r):.3f} p10 {np.percentile(r, 10):.3f} p90 {np.percentile(r, 90):.3f}"
            f" | n>1 only: mean {rm.mean():.3f} median {np.median(rm):.3f};  H/tol mean {(h / float(tol)).mean():.3f}; eq/tol mean {(eq / float(tol)).mean():.3f}"
        )

section("sanity: the Hausdorff distance can never exceed the equal-parameter distance (uniform-split methods)")
for name in ("fonttools_cu2qu", "kurbo_to_quads", "lyon", "harfbuzz", "reimpl_truong_presub"):
    for tol in ("1.0", "0.5", "0.1"):
        d = load(f"{name}.{tol}")
        h = np.maximum(d["haus_c2s"], d["haus_s2c"])
        print(f"  {name:22s} tol {tol}: max(H - eq) = {np.max(h - d['eq']):.3e} units; max eq/tol = {np.max(d['eq']) / float(tol):.9f}")

section("quadratics per cubic relative to fontTools cu2qu (same tolerance)")
for tol in ("1.0", "0.5", "0.1"):
    ref = load(f"fonttools_cu2qu.{tol}")["nseg"]
    for name in (
        "harfbuzz",
        "kurbo_to_quads",
        "lyon017",
        "cubic2quad_js",
        "skia",
        "skia_cal",
        "fontforge",
        "reimpl_truong_presub",
        "reimpl_colomitchi_half",
        "oracle_ramp",
        "oracle_truong",
        "oracle_qm",
        "oracle_tangent_implied",
        "oracle_tangent",
    ):
        p = os.path.join(B, "percubic", f"{name}.{tol}.txt")
        if not os.path.exists(p):
            continue
        n = load(f"{name}.{tol}")["nseg"]
        print(
            f"  tol {tol} {name:24s}: total {int(n.sum()):7d} ({n.sum() / ref.sum():.3f} x cu2qu); "
            f"fewer on {int((n < ref).sum()):6d}, same {int((n == ref).sum()):6d}, more on {int((n > ref).sum()):6d} cubics"
        )

section("extrema: cubics that start or end with a horizontal or vertical tangent at an extreme of their own bounding box")
# The node is an extremum of the cubic's box in the direction perpendicular to the tangent, and the cubic is
# monotone (control points inside the box of the end points), so its box is the box of its end points.
p0, p1, p2, p3 = cubics[:, 0], cubics[:, 1], cubics[:, 2], cubics[:, 3]
lo = np.minimum(p0, p3)
hi = np.maximum(p0, p3)
inside = ((p1 >= lo) & (p1 <= hi) & (p2 >= lo) & (p2 <= hi)).all(1)
hv_start = ((p1[:, 0] == p0[:, 0]) ^ (p1[:, 1] == p0[:, 1])) & ~zero_start
hv_end = ((p2[:, 0] == p3[:, 0]) ^ (p2[:, 1] == p3[:, 1])) & ~zero_end
sel = inside & (hv_start | hv_end) & ~straight
print(f"  {sel.sum()} cubics ({100 * sel.mean():.1f}%) are monotone with an H/V tangent at an end")
for name in ("fonttools_cu2qu", "harfbuzz", "kurbo_to_quads", "lyon", "cubic2quad_js", "skia", "skia_cal", "fontforge", "reimpl_truong_presub"):
    for tol in ("1.0", "0.5", "0.1"):
        p = os.path.join(B, "percubic", f"{name}.{tol}.txt")
        if not os.path.exists(p):
            continue
        d = load(f"{name}.{tol}")
        bb = d["bbox"][sel]
        print(
            f"  {name:22s} tol {tol}: spline box differs from the cubic's by > 0.001 units on {int((bb > 1e-3).sum()):6d} "
            f"({100 * (bb > 1e-3).mean():5.1f}%), > 0.1 on {int((bb > 0.1).sum()):6d}; mean {bb.mean():.4f}, max {bb.max():.3f}"
        )

section("how full is the tolerance budget: Hausdorff / tolerance")
for name in ("fonttools_cu2qu", "kurbo_to_quads", "cubic2quad_js", "skia", "skia_cal", "reimpl_truong_presub"):
    for tol in ("1.0", "0.5", "0.1"):
        d = load(f"{name}.{tol}")
        h = np.maximum(d["haus_c2s"], d["haus_s2c"]) / float(tol)
        print(
            f"  {name:22s} tol {tol}: mean {h.mean():.3f} median {np.median(h):.3f}; "
            f"share below 0.25: {100 * (h < 0.25).mean():5.1f}%  below 0.5: {100 * (h < 0.5).mean():5.1f}%  above 0.9: {100 * (h > 0.9).mean():4.1f}%"
        )
