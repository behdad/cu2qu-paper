#!/usr/bin/env python3
"""A fixed cu2qu workload, to be counted with `perf stat -e instructions:u`.

Usage: PYTHONPATH=<fontTools Lib> ft_perf.py cubics.json MODE
MODE: none    load only (the baseline to subtract)
      single  curve_to_quadratic on every cubic at tolerance 1.0, 3 times
      half    the same at tolerance 0.5
      multi   curves_to_quadratic on the cubics three at a time, 3 times
      random  cu2qu.benchmark's curves (0..2048, tolerance 0.05), single
"""
import json, random, sys
from fontTools.cu2qu.cu2qu import curve_to_quadratic, curves_to_quadratic

cubics = [[tuple(p) for p in c] for c in json.load(open(sys.argv[1]))]
mode = sys.argv[2]
groups = [cubics[i : i + 3] for i in range(0, len(cubics) - 2, 3)]
rnd = random.Random(1)
rand = [[tuple(float(rnd.randint(0, 2048)) for _ in range(2)) for _ in range(4)] for _ in range(20000)]
REPEAT = 3
if mode == "single":
    for _ in range(REPEAT):
        for c in cubics:
            curve_to_quadratic(c, 1.0)
elif mode == "half":
    for _ in range(REPEAT):
        for c in cubics:
            curve_to_quadratic(c, 0.5)
elif mode == "multi":
    tols = [1.0] * 3
    for _ in range(REPEAT):
        for g in groups:
            curves_to_quadratic(g, tols)
elif mode == "random":
    for _ in range(REPEAT):
        for c in rand:
            curve_to_quadratic(c, 0.05)
print(len(cubics) * REPEAT, len(groups) * REPEAT, len(rand) * REPEAT)
