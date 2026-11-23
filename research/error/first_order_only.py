#!/usr/bin/env python3
"""Is the first-order term alone ('error across the tangent', maximised exactly) a bound?
It is not: this measures by how much it falls short of the true distance on the corpus,
and how often a rule that used it alone would pick too few quadratics."""
import numpy as np
from common import *

c = load_cubics()
g = load_gt(W + '/gt_ramp.bin', 16)
H = np.fmax(g[:, :, 1], g[:, :, 2]); ep = g[:, :, 0]
fo_b, fo_s = g[:, :, 3], g[:, :, 4]
print(' n   cubics   across the spline tangent / H:  min    0.1%     1%     50%     99%  below H by >1%  |  across the cubic tangent / H:  min    0.1%     1%     50%     99%  below H by >1%')
for n in range(1, 7):
    m = np.isfinite(H[:, n - 1]) & (H[:, n - 1] > 1e-3)
    line = f'{n:2d}  {m.sum():6d}  '
    for fo in (fo_s, fo_b):
        r = fo[m, n - 1] / H[m, n - 1]
        q = np.percentile(r, [0.1, 1, 50, 99])
        line += f'                        {r.min():7.4f} ' + ' '.join(f'{x:7.4f}' for x in q) + f'    {100 * np.mean(r < 0.99):5.2f}%     |'
    print(line)
for tol in (1.0, 0.5, 0.1):
    nc = first_fit(H, tol * (1 + 1e-9))
    for name, fo in (('spline tangent', fo_s), ('cubic tangent', fo_b), ('larger of the two', np.fmax(fo_s, fo_b))):
        nf = first_fit(fo, tol)
        ok = nc > 0
        print(f'tol {tol}: rule by the first-order term across the {name:18s}: mean n {nf[ok].mean():.4f} (truth {nc[ok].mean():.4f}); too few quadratics for {(nf[ok] < nc[ok]).sum()} cubics')
