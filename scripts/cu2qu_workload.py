#!/usr/bin/env python3
"""A fixed cu2qu workload, for counting instructions and timing.

Usage: PYTHONPATH=<a fontTools Lib> cu2qu_workload.py cubics.json MODE [time]

MODE: none    load only (the baseline to subtract)
      single  curve_to_quadratic on every cubic at tolerance 1.0
      half    the same at tolerance 0.5
      multi   curves_to_quadratic on the cubics three at a time, tolerance 1.0

Each mode goes through the corpus REPEAT times. With "time" as a third
argument it prints the best time per call over 7 runs instead (microseconds).
The last line printed is: mode, calls, compiled.
"""
import json
import sys
import time

from fontTools.cu2qu import cu2qu
from fontTools.cu2qu.cu2qu import curve_to_quadratic, curves_to_quadratic

REPEAT = 3
cubics = [[tuple(float(v) for v in p) for p in c] for c in json.load(open(sys.argv[1]))]
mode = sys.argv[2]
timing = len(sys.argv) > 3 and sys.argv[3] == "time"
groups = [cubics[i : i + 3] for i in range(0, len(cubics) - 2, 3)]
tolerances = [1.0] * 3


def run():
    if mode == "single":
        for c in cubics:
            curve_to_quadratic(c, 1.0)
        return len(cubics)
    if mode == "half":
        for c in cubics:
            curve_to_quadratic(c, 0.5)
        return len(cubics)
    if mode == "multi":
        for g in groups:
            curves_to_quadratic(g, tolerances)
        return len(groups)
    return 0


if timing:
    best = float("inf")
    calls = 1
    for _ in range(7):
        start = time.perf_counter()
        calls = run()
        best = min(best, time.perf_counter() - start)
    print("%.4f" % (best / max(calls, 1) * 1e6))
else:
    calls = 0
    for _ in range(REPEAT):
        calls += run()
print(mode, calls, cu2qu.COMPILED)
