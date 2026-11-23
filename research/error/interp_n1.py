#!/usr/bin/env python3
"""The interpolation caveat for n = 1.  Pairs of masters from interp_experiment.py for which
BOTH masters are within tol of their own tangent-intersection quadratic by the true distance.
The instance half-way has the cubic (A + B) / 2 and the quadratic with control point
(q_A + q_B) / 2.  Its Hausdorff distance, by brute force on dense samples (all pairs;
accuracy about 1e-3 of the curve's sampling step), against tol."""
import numpy as np
from common import *
from closed import n1_setup

out = W + '/interp'
NMAX = 10
A = load_cubics(f'{out}/A.bin'); B = load_cubics(f'{out}/B.bin'); M = load_cubics(f'{out}/M.bin')
gA = load_gt(f'{out}/A.gt', NMAX); gB = load_gt(f'{out}/B.gt', NMAX)
HA = np.fmax(gA[:, 0, 1], gA[:, 0, 2]); HB = np.fmax(gB[:, 0, 1], gB[:, 0, 2])
eA = gA[:, 0, 0]; eB = gB[:, 0, 0]
qA = n1_setup(A)[0]; qB = n1_setup(B)[0]

def hausdorff(cub, q, T=1501):
    t = np.linspace(0, 1, T)[:, None]
    Bc = (1-t)**3*cub[0] + 3*(1-t)**2*t*cub[1] + 3*(1-t)*t*t*cub[2] + t**3*cub[3]
    Q = (1-t)**2*cub[0] + 2*(1-t)*t*q + t*t*cub[3]
    D = np.sqrt(((Bc[:, None, :] - Q[None, :, :]) ** 2).sum(-1))
    return max(D.min(0).max(), D.min(1).max())

print('tol   pairs with both masters a single quadratic: by equal-parameter   by true distance   instance outside tol (true-distance pairs)   median / 90% / max of H_mid/tol among them')
for tol in (1.0, 0.5):
    both_ep = (eA <= tol) & (eB <= tol)
    both = (HA <= tol) & (HB <= tol) & np.isfinite(qA).all(1) & np.isfinite(qB).all(1)
    idx = np.where(both)[0][:3000]
    r = np.array([hausdorff(M[i], (qA[i] + qB[i]) / 2) for i in idx]) / tol
    viol = r > 1.002
    # the same for the pairs today's test accepts: must never violate
    idx2 = np.where(both_ep)[0][:1000]
    r2 = np.array([hausdorff(M[i], (qA[i] + qB[i]) / 2) for i in idx2]) / tol if len(idx2) else np.array([0.0])
    print(f'{tol:4.1f}   {both_ep.sum():6d}                                                {both.sum():6d}            {viol.sum():5d} of {len(idx)} ({100 * viol.mean():.1f}%)'
          f'                       {np.median(r):.2f} / {np.percentile(r, 90):.2f} / {r.max():.2f}     (equal-parameter pairs: max {r2.max():.3f})')
