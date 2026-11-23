#!/usr/bin/env python3
"""The tables of the report, random and adversarial sets (20,000 cubics each, n = 1..8).

Table 4: for every set, over all (cubic, n) pairs: how often the closed form applies,
the number of under-estimates (bound < measured Hausdorff distance, beyond 1e-6
relative + 1e-9 of the coordinate size), and the ratios (a)/(c), (b)/(c).
Table 5: the quadratics each measure picks.  Tolerances are scaled to the size of
each cubic's control polygon so that the wild sets are tested where n <= 8 matters:
tol = size / 100 and size / 400 (cubics that the equal-parameter test cannot fit
with 8 quadratics are left out of the means)."""
import sys
import numpy as np
from common import *
from closed import *

NMAX = 8
SETS = sys.argv[1] if len(sys.argv) > 1 else 'sets'       # subdirectory of W
names = ['random', 'fontlike', 'arcs', 'cusp', 'loop', 'inflect', 'coincident', 'straight', 'parallel', 'tiny_huge']
tot_under = 0; tot_pairs = 0; tot_valid = 0
rows4 = []; rows5 = []
for name in names:
    c = load_cubics(f'{W}/{SETS}/{name}.bin')
    g = load_gt(f'{W}/{SETS}/{name}.gt_ramp', NMAX)
    H = np.fmax(g[:, :, 1], g[:, :, 2]); ep = g[:, :, 0]
    N = len(c)
    scale = np.abs(c).max((1, 2))
    size = np.maximum(np.ptp(c[:, :, 0], axis=1), np.ptp(c[:, :, 1], axis=1))
    B = np.full((N, NMAX), np.inf)
    under = 0; valid_pairs = 0; pairs = 0; raw = []; used = []; ra = []
    worst = None
    for n in range(1, NMAX + 1):
        b, v = n1_closed(c) if n == 1 else ramp_closed(c, n)
        h = H[:, n - 1]; e = ep[:, n - 1]
        m = np.isfinite(h) & np.isfinite(e)
        B[:, n - 1] = np.where(v, np.minimum(b, e), e)
        u = m & v & (b < h * (1 - 1e-6) - 1e-9 * scale)
        if u.any() and worst is None:
            i = np.where(u)[0][0]
            worst = (n, i, c[i].tolist(), b[i], h[i])
        under += u.sum()
        valid_pairs += (m & v).sum(); pairs += m.sum()
        mm = m & (h > 1e-7 * scale)
        mv = mm & v
        raw.append(b[mv] / h[mv]); used.append(B[mm, n - 1] / h[mm]); ra.append(e[mm] / h[mm])
    tot_under += under; tot_valid += valid_pairs; tot_pairs += pairs
    raw = np.concatenate(raw) if valid_pairs else np.array([np.nan]); used = np.concatenate(used); ra = np.concatenate(ra)
    rows4.append(f'{name:11s} {pairs:7d}  {100 * valid_pairs / pairs:5.1f}%  {under:6d}   |  {np.median(ra):7.3f} {np.percentile(ra, 90):8.2f}   |  '
                 f'{np.nanmin(raw):9.6f} {np.nanmedian(raw):8.4f} {np.nanpercentile(raw, 90):8.3f} {np.nanpercentile(raw, 99):9.2f}   |  {np.median(used):7.3f} {np.percentile(used, 90):8.2f}')
    if worst:
        rows4.append(f'      first under-estimate: n {worst[0]} cubic {worst[1]} {worst[2]} bound {worst[3]:.8g} measured {worst[4]:.8g}')
    line = f'{name:11s}'
    for div in (100.0, 400.0):
        tol = size / div
        na = first_fit(ep / tol[:, None], 1.0); nb = first_fit(B / tol[:, None], 1.0); nc = first_fit(H / tol[:, None], 1 + 1e-9)
        ok = (na > 0) & (size > 0)
        line += f'   {ok.sum():6d}  {na[ok].mean():6.3f} {nb[ok].mean():6.3f} {nc[ok].mean():6.3f}   {(nb[ok] < nc[ok]).sum():4d}     |'
    rows5.append(line)

print('Table 4.  Random and adversarial sets, ramp, n = 1..8, all (cubic, n) pairs.')
print('set           pairs   valid   under   |  (a)/(c): 50%     90%   |  (b)/(c) where the closed form is valid: min  50%  90%  99%   |  (b)/(c) as used (with fallback): 50%  90%')
print('\n'.join(rows4))
print(f'total: {tot_pairs} pairs measured, closed form valid for {tot_valid}, under-estimates {tot_under}')
print()
print('Table 5.  The same sets: mean quadratics per cubic, tolerance relative to the size of the control polygon.')
print('set          tol = size/100: cubics   (a)    (b)    (c)  (b)<(c)  |  tol = size/400: cubics   (a)    (b)    (c)  (b)<(c)')
print('\n'.join(rows5))
