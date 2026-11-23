#!/usr/bin/env python3
"""n = 1 rule of thumb: the distance between a cubic and its tangent-intersection quadratic
is about (3/8) |l1 + l2 - 4/3| h, with h the distance of the tangent intersection from the
chord and l1, l2 the handle lengths as fractions of the tangent legs.  (The leading term of
the bound, at t = 1/2, for legs of similar length.)  How good is it on the corpus?"""
import numpy as np
from common import *
from closed import n1_setup, cross, dot

c = load_cubics()
g = load_gt(W + '/gt_ramp.bin', 16)
H1 = np.fmax(g[:, 0, 1], g[:, 0, 2])
q, a, z, l1, l2 = n1_setup(c)
chord = c[:, 3] - c[:, 0]
with np.errstate(divide='ignore', invalid='ignore'):
    h = np.abs(cross(a, z)) / np.sqrt(dot(chord, chord))
    est = 3 / 8 * np.abs(l1 + l2 - 4 / 3) * h
m = np.isfinite(est) & np.isfinite(H1) & (l1 > 0) & (l1 < 2) & (l2 > 0) & (l2 < 2) & (H1 > 1e-3)
r = est[m] / H1[m]
print('cubics:', m.sum())
print('estimate / measured distance: percentiles 1/5/25/50/75/95/99:', ' '.join(f'{x:.3f}' for x in np.percentile(r, [1, 5, 25, 50, 75, 95, 99])))
for tol in (1.0, 0.5):
    sel = m & (H1 < 3 * tol) & (H1 > tol / 3)
    rr = est[sel] / H1[sel]
    print(f'  where the distance is within a factor 3 of tol = {tol}: {sel.sum()} cubics, percentiles 5/25/50/75/95: ' + ' '.join(f'{x:.3f}' for x in np.percentile(rr, [5, 25, 50, 75, 95])))
