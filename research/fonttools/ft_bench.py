#!/usr/bin/env python3
"""Time fontTools' cu2qu on cubics from fonts and on its own benchmark's
random curves.

Usage: PYTHONPATH=<fontTools Lib> ft_bench.py cubics.json
Prints microseconds per call (best of several repeats).
"""
import json, random, sys, time
from fontTools.cu2qu import cu2qu
from fontTools.cu2qu.cu2qu import curve_to_quadratic, curves_to_quadratic

cubics = [[tuple(p) for p in c] for c in json.load(open(sys.argv[1]))]


def best(f, repeats=7):
    t = float("inf")
    for _ in range(repeats):
        t0 = time.perf_counter()
        f()
        t = min(t, time.perf_counter() - t0)
    return t


def font_single(tol):
    def run():
        for c in cubics:
            curve_to_quadratic(c, tol)
    return best(run) / len(cubics) * 1e6


groups = [cubics[i : i + 3] for i in range(0, len(cubics) - 2, 3)]


def font_multi(tol):
    tols = [tol] * 3
    def run():
        for g in groups:
            curves_to_quadratic(g, tols)
    return best(run) / len(groups) * 1e6


# The curves of fontTools/cu2qu/benchmark.py: random points in 0..2048,
# tolerance 0.05.
rnd = random.Random(1)
rand = [[tuple(float(rnd.randint(0, 2048)) for _ in range(2)) for _ in range(4)] for _ in range(3000)]
rgroups = [rand[i : i + 3] for i in range(0, len(rand), 3)]


def rand_single():
    def run():
        for c in rand:
            curve_to_quadratic(c, 0.05)
    return best(run) / len(rand) * 1e6


def rand_multi():
    tols = [0.05] * 3
    def run():
        for g in rgroups:
            curves_to_quadratic(g, tols)
    return best(run) / len(rgroups) * 1e6


n_quads = sum(len(curve_to_quadratic(c, 1.0)) - 2 for c in cubics)
print("compiled" if cu2qu.COMPILED else "pure", end=" ")
print(
    "font_single_1.0 %.3f font_single_0.5 %.3f font_multi3_1.0 %.3f random_single_0.05 %.3f random_multi3_0.05 %.3f quads_per_cubic_1.0 %.4f"
    % (font_single(1.0), font_single(0.5), font_multi(1.0), rand_single(), rand_multi(), n_quads / len(cubics))
)
