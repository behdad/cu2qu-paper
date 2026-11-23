#!/usr/bin/env python3
"""First look: how much is there to gain?  n by the equal-parameter error (a),
by the true Hausdorff distance (c), and by the first-order 'error across the
tangent' measured by dense sampling (not a bound; the target of the formulas)."""
import numpy as np
from common import *

c = load_cubics()
g = load_gt(W + '/gt_ramp.bin', 16)
ep, hsb, hbs, fob, fos = (g[:, :, j] for j in range(5))
H = np.fmax(hsb, hbs)
print('fraction of (cubic, n) where h_sb > h_bs:', np.nanmean(hsb > hbs))
for tol in (1.0, 0.5, 0.1):
    n_a = first_fit(ep, tol)
    n_c = first_fit(H, tol)
    n_fb = first_fit(fob, tol)
    n_fs = first_fit(fos, tol)
    n_fm = first_fit(np.fmax(fob, fos), tol)
    ok = (n_a > 0)
    print(f'tol {tol}: N={ok.sum()}  mean n: equal-param {n_a[ok].mean():.4f}  true {n_c[ok].mean():.4f} ({100*(1-n_c[ok].mean()/n_a[ok].mean()):.1f}% fewer)'
          f'  first-order(B) {n_fb[ok].mean():.4f}  first-order(S) {n_fs[ok].mean():.4f}  first-order(max) {n_fm[ok].mean():.4f}')
    # joint histogram
    mx = int(n_a[ok].max())
    print('   n_a \\ n_true ' + ' '.join(f'{j:6d}' for j in range(1, 9)))
    for a in range(1, min(mx, 10) + 1):
        row = [(ok & (n_a == a) & (n_c == j)).sum() for j in range(1, 9)]
        print(f'   {a:2d}            ' + ' '.join(f'{v:6d}' for v in row))
    # first-order measure vs truth: where does it under-call n?
    print('   first-order(max) n < true n:', (n_fm[ok] < n_c[ok]).sum(), ' > true n:', (n_fm[ok] > n_c[ok]).sum())

# ratios at the n chosen by equal-parameter test, tol 1
for tol in (1.0, 0.5):
    n_a = first_fit(ep, tol)
    ok = n_a >= 1
    idx = np.arange(len(c))[ok]
    e = ep[idx, n_a[ok] - 1]; h = H[idx, n_a[ok] - 1]; f = np.fmax(fob, fos)[idx, n_a[ok] - 1]
    m = h > 1e-9
    r = e[m] / h[m]
    q = np.percentile(r, [0, 5, 25, 50, 75, 95, 100])
    print(f'tol {tol}: ep/H at the n picked by ep: percentiles 0/5/25/50/75/95/100:', ' '.join(f'{x:.3f}' for x in q))
    r = f[m] / h[m]
    q = np.percentile(r, [0, 0.1, 1, 5, 25, 50, 75, 95, 99, 99.9, 100])
    print(f'          first-order/H:  percentiles 0/.1/1/5/25/50/75/95/99/99.9/100:', ' '.join(f'{x:.4f}' for x in q))
