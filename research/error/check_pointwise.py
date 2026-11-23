#!/usr/bin/env python3
"""n >= 2, ramp, corpus: the pointwise formulas maximised by dense sampling (U points per
quadratic) against the measured Hausdorff distance.

  first = |phi_k(u)| |d x V(u)| / |V(u)|          the error across the spline's tangent
  exact = |E.n - (C/|V|) eps^2|                   the distance along the normal line (Lemma N)

Is the maximum of `exact` the Hausdorff distance?  (Cubics where Lemma N does not apply
somewhere, y < -1, are counted separately.)   Usage: check_pointwise.py [U]"""
import sys
import numpy as np
from common import *
from geo import *

U = int(sys.argv[1]) if len(sys.argv) > 1 else 129
c = load_cubics()
NMAX = 16
g = load_gt(W + '/gt_ramp.bin', NMAX)
H = np.fmax(g[:, :, 1], g[:, :, 2])
print(f'{U} points per quadratic')
print(' n   cubics  y<-1 |  first/H:  min     0.1%     1%     50%     99%    max  |  exact/H:  min       0.1%      1%      50%      99%     max')
for n in range(2, 9):
    m = np.isfinite(H[:, n - 1]) & (H[:, n - 1] > 1e-6)
    if not m.any():
        continue
    first, exact, bad = pointwise_dense(c[m], weights('ramp', n), U=U)
    h = H[m, n - 1][~bad]
    r1 = first[~bad] / h; r2 = exact[~bad] / h
    print(f'{n:2d}  {m.sum():6d}  {bad.sum():3d}  |          {r1.min():.4f} ' + ' '.join(f'{x:.4f}' for x in np.percentile(r1, [0.1, 1, 50, 99])) + f' {r1.max():6.3f}  |'
          f'          {r2.min():.6f} ' + ' '.join(f'{x:.6f}' for x in np.percentile(r2, [0.1, 1, 50, 99])) + f' {r2.max():7.4f}')
