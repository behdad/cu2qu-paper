#!/usr/bin/env python3
"""Q4 (and Q2, Q3): canonical shapes, n = 1..6, geometric (Hausdorff) error of
each closed-form scheme against numerically optimised splines.

Columns (all Hausdorff distance in font units unless said otherwise):
  ep-bound  Phi(n) |d3| / n^3 for the cu2qu ramp (what cu2qu tests against)
  ramp      cu2qu's spline (n = 1: tangent intersection)
  optu      d3 family, equal-parameter-optimal weights with tangents kept
  T-uni     tangent-intersection spline at equal parameter knots
  T-aff     tangent-intersection spline at equal affine-arc-length knots
  w*        best weights of the d3 family (tangents kept) for the geometric distance
  ws*       w* plus sliding the two end points along the end tangents
  TAN*      best implied-on-curve spline with both end tangents kept (yardstick)
  zero      d3 family with s = 0 (tangents not kept)
  FREE*     best implied-on-curve spline, tangents not kept
"""
import sys
import numpy as np
import chd, schemes
from opt import best_spline
from common import phi_max, ramp as ramp_w, power_basis

SHAPES = {
    'quarter circle R=500':       [(500, 0), (500, 276.142), (276.142, 500), (0, 500)],
    'corpus bowl arc':            [(300, 0), (447, 0), (547, 63), (547, 192)],
    'superellipse-like quarter':  [(500, 0), (500, 330), (330, 500), (0, 500)],
    'S-curve':                    [(0, 0), (200, 0), (100, 300), (300, 300)],
    'flat arc, long handles':     [(0, 0), (450, 30), (550, 30), (1000, 0)],
    'flat, one-sided handles':    [(0, 0), (600, 0), (800, 20), (1000, 60)],
    'near-cusp':                  [(0, 0), (270, 300), (30, 300), (300, 0)],
    'hook (line into curve)':     [(0, 0), (0, 300), (200, 400), (500, 400)],
}


def row(c, n):
    c = np.asarray(c, dtype=float)
    C = c[None]
    b3 = power_basis(c)[3]
    d3 = float(np.hypot(*b3))
    out = {}
    out['ep-bound'] = (phi_max(ramp_w(n)) if n > 1 else float('nan')) * d3 / n ** 3
    Tu = schemes.tangents(C, schemes.uniform_knots(1, n), return_valid=True)
    Ta = schemes.tangents(C, schemes.affine_knots(C, n), return_valid=True)
    Qr = schemes.ramp(C, n)[0] if n > 1 else Tu[0][0]
    Qo = schemes.optu(C, n)[0] if n > 1 else Tu[0][0]
    Qz = schemes.zero(C, n)[0]
    h = lambda Q: chd.hd1(c, Q, m=64)[0]
    out['ramp'] = h(Qr) if (n > 1 or Tu[1].all()) else float('nan')
    out['optu'] = h(Qo) if (n > 1 or Tu[1].all()) else float('nan')
    out['T-uni'] = h(Tu[0][0]) if Tu[1].all() else float('nan')
    out['T-aff'] = h(Ta[0][0]) if Ta[1].all() else float('nan')
    out['zero'] = h(Qz)
    starts_t = [Qr, Qo] + ([Tu[0][0]] if Tu[1].all() else []) + ([Ta[0][0]] if Ta[1].all() else [])
    for mode, key, starts in (('w', 'w*', [Qo, Qr]), ('ws', 'ws*', [Qo, Qr]),
                              ('tan', 'TAN*', starts_t), ('free', 'FREE*', starts_t + [Qz])):
        if n == 1 and mode != 'free':
            out[key] = out['T-uni']
            continue
        r = best_spline(c, n, mode, starts)
        out[key] = r[0] if r is not None else float('nan')
    return out


if __name__ == '__main__':
    cols = ['ep-bound', 'ramp', 'optu', 'T-uni', 'T-aff', 'w*', 'ws*', 'TAN*', 'zero', 'FREE*']
    names = sys.argv[1:] or list(SHAPES)
    for name in names:
        c = SHAPES[name]
        b3 = power_basis(np.asarray(c, dtype=float))[3]
        print(f'== {name}: {c}   |d3| = {np.hypot(*b3):.1f}')
        print('   n  ' + ''.join(f'{k:>10s}' for k in cols))
        for n in range(1, 7):
            r = row(c, n)
            print(f'  {n:2d}  ' + ''.join(f'{r[k]:10.4f}' for k in cols), flush=True)
        print()
