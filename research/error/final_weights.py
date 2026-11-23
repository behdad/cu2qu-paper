#!/usr/bin/env python3
"""The other weightings, corpus: s = 0 (every off-curve point at the 'midpoint' control;
end tangents not kept) and the alternating +1,-1,... (even n only; n = 1 is the
tangent intersection).  Mean quadratics per cubic by (a) the equal-parameter error,
(b) the closed-form geometric bound (spline_closed), (c) the measured Hausdorff distance."""
import numpy as np
from common import *
from closed import *
from geo import weights

c = load_cubics()
NMAX = 16
N = len(c)
EPS = 1e-9
scale = np.abs(c).max((1, 2))
res = {}
for mode, fn in (('ramp', 'gt_ramp.bin'), ('zero', 'gt_zero.bin'), ('alt', 'gt_alt.bin')):
    g = load_gt(f'{W}/{fn}', NMAX)
    H = np.fmax(g[:, :, 1], g[:, :, 2]); ep = g[:, :, 0]
    B = np.full((N, NMAX), np.inf)
    under = 0; pairs = 0
    for n in range(1, NMAX + 1):
        h = H[:, n - 1]; e = ep[:, n - 1]
        if not np.isfinite(h).any():
            continue
        if n == 1 and mode != 'zero':
            b, v = n1_closed(c)
        elif mode == 'ramp':
            b, v = ramp_closed(c, n)
        elif mode == 'alt' and n % 2:
            continue
        else:
            b, v = spline_closed(c, n, weights(mode, n) if n > 1 else np.zeros(1))
        m = np.isfinite(h) & np.isfinite(e)
        B[:, n - 1] = np.where(v, np.minimum(b, e), e)
        B[~np.isfinite(e), n - 1] = np.inf
        under += (m & v & (b < h * (1 - 1e-6) - 1e-9 * scale)).sum(); pairs += (m & v).sum()
    res[mode] = (ep, B, H)
    print(f'{mode}: closed form valid for {pairs} (cubic, n) pairs, under-estimates {under}')

print()
print('Mean quadratics per cubic, corpus.  Savings are against the ramp with the equal-parameter test (cu2qu today).')
print('weighting                         measure                  tol 1.0           tol 0.5           tol 0.1')
base = {}
for mode, label in (('ramp', 'ramp (cu2qu)'), ('alt', 'alternating +1,-1 (even n)'), ('zero', 's = 0 (tangents not kept)')):
    ep, B, H = res[mode]
    for name, arr in (('(a) equal-parameter', ep), ('(b) closed form', B), ('(c) true least n', H)):
        line = f'{label:33s} {name:22s}'
        for tol in (1.0, 0.5, 0.1):
            nn = first_fit(arr, tol * (1 + EPS) if arr is H else tol)
            ok = first_fit(res['ramp'][0], tol) > 0
            assert (nn[ok] > 0).all(), (mode, name, tol, (nn[ok] == 0).sum())
            mean = nn[ok].mean()
            base.setdefault(tol, mean)
            line += f'   {mean:.4f} ({100 * (1 - mean / base[tol]):5.1f}%)'
            if arr is B:
                nc = first_fit(H, tol * (1 + EPS))
                assert (nn[ok] >= nc[ok]).all()
        print(line)
