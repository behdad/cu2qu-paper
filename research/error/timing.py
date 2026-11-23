#!/usr/bin/env python3
"""Rough cost of the rule, in the same language as the reference: pure-Python fontTools
curve_to_quadratic (equal-parameter test by trial) against cu2qu_geo.quadratics_needed
(equal-parameter error in closed form, then the geometric bound), on 5,000 corpus cubics.
The second builds no output spline, so add the cost of building one spline of n quadratics."""
import sys, time
sys.path.insert(0, '/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad/ft-old/Lib')
import numpy as np
from common import *
import cu2qu_geo as G
from fontTools.cu2qu import cu2qu as ft
from fontTools.cu2qu.cu2qu import curve_to_quadratic

assert not ft.COMPILED, 'expected the pure-Python fontTools'
c = load_cubics()
idx = np.random.default_rng(3).choice(len(c), 5000, replace=False)
tup = [[tuple(p) for p in c[i]] for i in idx]
cpx = [[complex(*p) for p in c[i]] for i in idx]
for tol in (1.0, 0.5, 0.1):
    best_ft = best_geo = best_build = 1e9
    for rep in range(3):
        t0 = time.perf_counter(); n_ft = [len(curve_to_quadratic(x, tol)) - 2 for x in tup]; t1 = time.perf_counter()
        n_geo = [G.quadratics_needed(x, tol) for x in cpx]; t2 = time.perf_counter()
        for x, n in zip(cpx, n_geo):
            if n > 1:
                G.ramp_spline(*x, n)
        t3 = time.perf_counter()
        best_ft = min(best_ft, t1 - t0); best_geo = min(best_geo, t2 - t1); best_build = min(best_build, t3 - t2)
    print(f'tol {tol}: fontTools (pure Python) {1e6 * best_ft / len(idx):6.1f} us/cubic, mean n {np.mean(n_ft):.3f};   '
          f'geometric rule {1e6 * best_geo / len(idx):6.1f} us/cubic (+ {1e6 * best_build / len(idx):4.1f} to build the spline), mean n {np.mean(n_geo):.3f}')
