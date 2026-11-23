#!/usr/bin/env python3
"""Compare two runs of the measuring instrument (e.g. scan-based vs exact nearest point,
or two resolutions): relative differences of the Hausdorff distance."""
import sys
import numpy as np
from common import *

a = load_gt(sys.argv[1], int(sys.argv[3]))
b = load_gt(sys.argv[2], int(sys.argv[3]))
Ha = np.fmax(a[:, :, 1], a[:, :, 2]); Hb = np.fmax(b[:, :, 1], b[:, :, 2])
m = np.isfinite(Ha) & np.isfinite(Hb)
print('entries compared:', m.sum(), ' only in first:', (np.isfinite(Ha) & ~np.isfinite(Hb)).sum(), ' only in second:', (~np.isfinite(Ha) & np.isfinite(Hb)).sum())
d = np.abs(Ha[m] - Hb[m])
rel = d / np.maximum(np.maximum(Ha[m], Hb[m]), 1e-12)
big = Ha[m] > 1e-3
print('abs diff: max %.3g  99.9%% %.3g  median %.3g' % (d.max(), np.percentile(d, 99.9), np.median(d)))
print('rel diff where H > 1e-3: max %.3g  99.99%% %.3g  99.9%% %.3g  99%% %.3g  median %.3g' % (rel[big].max(), np.percentile(rel[big], 99.99), np.percentile(rel[big], 99.9), np.percentile(rel[big], 99), np.median(rel[big])))
print('first > second by more than 1e-6 rel:', (big & (Ha[m] > Hb[m] * (1 + 1e-6))).sum(), '  second > first:', (big & (Hb[m] > Ha[m] * (1 + 1e-6))).sum())
ii = np.argwhere(m)
o = np.argsort(-d)[:5]
for j in o:
    i, n = ii[j]
    print('  cubic', i, 'n', n + 1, 'H first', Ha[i, n], 'second', Hb[i, n])
