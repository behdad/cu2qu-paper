#!/usr/bin/env python3
"""A caveat for interpolatable masters.  The equal-parameter error is convex in the
control points, so if two masters fit with the same n, every interpolated instance
fits.  The geometric distance is not.  Experiment: pairs of masters (a corpus cubic
and a perturbed copy), n = the least n >= 2 for which both masters are within tol by
the *true* geometric distance, then the instance half-way: how often is it outside
tol, and by how much?  (n >= 2 only: there the spline of the instance is the
interpolation of the masters' splines.)

Step 1 (no argument) writes the three sets; run gt on them; step 2 (argument 'report') reads the results."""
import sys, subprocess
import numpy as np
from common import *

out = W + '/interp'
NMAX = 10
if len(sys.argv) < 2:
    import os
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(7)
    c = load_cubics()
    A = c[rng.choice(len(c), 20000, replace=False)]
    chord = np.hypot(*(A[:, 3] - A[:, 0]).T)[:, None, None]
    B = A + rng.normal(size=A.shape) * 0.08 * chord        # a different master: every point moved by ~8% of the chord
    M = (A + B) / 2
    for name, x in (('A', A), ('B', B), ('M', M)):
        x.reshape(-1, 8).tofile(f'{out}/{name}.bin')
        subprocess.run([W + '/gt', f'{out}/{name}.bin', f'{out}/{name}.gt', str(NMAX), '0', '0', '64', '256'], check=True)
    sys.exit(0)

gA, gB, gM = (load_gt(f'{out}/{n}.gt', NMAX) for n in 'ABM')
HA, HB, HM = (np.fmax(g[:, :, 1], g[:, :, 2]) for g in (gA, gB, gM))
eA, eB, eM = (g[:, :, 0] for g in (gA, gB, gM))
print('tol   pairs   n by equal-parameter (both masters)   n by true distance (both masters)   instance outside tol: equal-parameter   true-distance rule   worst H_mid/tol   90th pct of violators')
for tol in (1.0, 0.5):
    def least(XA, XB):
        ok = (np.fmax(XA, XB)[:, 1:] <= tol * (1 + 1e-9))
        n = np.argmax(ok, axis=1) + 2
        n[~ok.any(1)] = 0
        return n
    n_ep = least(eA, eB); n_geo = least(HA, HB)
    m = (n_ep > 0) & (n_geo > 0)
    i = np.arange(len(m))[m]
    h_mid_ep = HM[i, n_ep[m] - 1]; h_mid_geo = HM[i, n_geo[m] - 1]
    viol_ep = h_mid_ep > tol * (1 + 1e-9); viol = h_mid_geo > tol * (1 + 1e-9)
    print(f'{tol:4.1f}  {m.sum():6d}   {n_ep[m].mean():.3f}                                 {n_geo[m].mean():.3f}                              '
          f'{viol_ep.sum():6d}                      {viol.sum():6d} ({100 * viol.mean():.1f}%)       {(h_mid_geo / tol).max():.2f}            {np.percentile(h_mid_geo[viol] / tol, 90) if viol.any() else 0:.2f}')
