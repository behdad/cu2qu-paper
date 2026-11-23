#!/usr/bin/env python3
"""Q2: what do the numerically best tangent-kept splines look like?

For a few shapes and n = 2..5: the parameter of the cubic point nearest each
implied on-curve point, against the knots of equal affine arc length; the
signed offset of the on-curve points; the largest and smallest signed error of
each quadratic; and the error of the plain tangent-intersection spline at
those same knots (it should be about twice the optimum: the optimum is that
spline moved inward by half its one-sided error)."""
import sys
import numpy as np
import chd, schemes
from opt import best_spline
from canon import SHAPES
from common import power_basis

def foot(c, P):
    # nearest point on cubic to P: dense search + Newton via chd.resid on a fake 1-quad spline is awkward; do numpy
    ts = np.linspace(0, 1, 20001)
    b0, b1, b2, b3 = power_basis(c)
    B = b0 + ts[:, None] * (b1 + ts[:, None] * (b2 + ts[:, None] * b3))
    D = b1 + ts[:, None] * (2 * b2 + 3 * ts[:, None] * b3)
    j = ((B - P) ** 2).sum(-1).argmin()
    v = P - B[j]
    nrm = np.array([-D[j, 1], D[j, 0]]) / np.hypot(*D[j])
    return ts[j], float(v @ nrm)

for name in sys.argv[1:] or ['quarter circle R=500', 'corpus bowl arc', 'hook (line into curve)', 'S-curve']:
    c = np.asarray(SHAPES[name], dtype=float)
    C = c[None]
    print('==', name)
    for n in (2, 3, 4, 5):
        Tu = schemes.tangents(C, schemes.uniform_knots(1, n), return_valid=True)
        ka = schemes.affine_knots(C, n)
        Ta = schemes.tangents(C, ka, return_valid=True)
        starts = [schemes.ramp(C, n)[0], schemes.optu(C, n)[0]] + ([Tu[0][0]] if Tu[1].all() else []) + ([Ta[0][0]] if Ta[1].all() else [])
        h, Q, si = best_spline(c, n, 'tan', starts)
        on = [c[0]] + [(Q[k - 1] + Q[k]) / 2 for k in range(1, n)] + [c[3]]
        feet = [foot(c, P) for P in on]
        r, gx, gy, tf = chd.resid(c, Q, 32)
        # extrema of r per quad
        ext = []
        for k in range(n):
            seg = r[k * 32:(k + 1) * 32 + 1]
            ext.append((seg.max(), seg.min()))
        print(f'  n = {n}: TAN* = {h:.4f} (start {si})')
        print('     feet t of on-curve points : ' + ' '.join(f'{t:.4f}' for t, _ in feet))
        print('     equal affine-arc knots    : ' + ' '.join(f'{t:.4f}' for t in ka[0]))
        print('     signed offset of on-curve : ' + ' '.join(f'{d:+.4f}' for _, d in feet))
        print('     per-quad max / min of signed error: ' + '  '.join(f'[{a:+.4f},{b:+.4f}]' for a, b in ext))
        # tangent-intersection spline at the feet knots
        kn = np.array([[t for t, _ in feet]])
        Tf = schemes.tangents(C, kn, return_valid=True)
        if Tf[1].all():
            print(f'     T at those knots: {chd.hd1(c, Tf[0][0], m=48)[0]:.4f}   T at affine knots: {chd.hd1(c, Ta[0][0], m=48)[0]:.4f}')
