#!/usr/bin/env python3
"""Sanity checks of the measuring instrument gt.c on the corpus:
 1. its equal-parameter error for n >= 2 against the closed form |d3| Phi(n) / n^3;
 2. the n fontTools picks against the first n whose equal-parameter error fits."""
import sys
sys.path.insert(0, '/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad/ft-old/Lib')
import numpy as np
from common import *
from fontTools.cu2qu.cu2qu import curve_to_quadratic

c = load_cubics()
g = load_gt(W + '/gt_ramp.bin', 16)
d3 = np.hypot(*d3_of(c).T)
worst = 0
for n in range(2, 17):
    ep = g[:, n - 1, 0]
    m = np.isfinite(ep)
    pred = d3[m] * Phi_ramp(n) / n ** 3
    rel = np.abs(ep[m] - pred) / np.maximum(pred, 1e-9)
    worst = max(worst, rel.max()) if m.any() else worst
    print(f'n={n:2d}  computed {m.sum():6d}  max rel diff of ep vs closed form {rel.max() if m.any() else 0:.2e}')
print('worst', worst)

for tol in (1.0, 0.5, 0.1):
    n_ep = first_fit(g[:, :, 0], tol)
    bad = 0
    hist = {}
    for i in range(len(c)):
        try:
            sp = curve_to_quadratic([tuple(p) for p in c[i]], tol)
            n_ft = len(sp) - 2
        except Exception:
            n_ft = 0
        hist[n_ft] = hist.get(n_ft, 0) + 1
        if n_ft != n_ep[i]:
            bad += 1
            if bad <= 5:
                print('  mismatch', i, c[i].tolist(), 'fontTools', n_ft, 'gt', n_ep[i], g[i, :4, 0])
    print(f'tol {tol}: fontTools n vs first n with ep <= tol: {bad} mismatches of {len(c)};  n histogram {dict(sorted(hist.items()))};  mean {sum(k*v for k,v in hist.items())/len(c):.4f}')
