#!/usr/bin/env python3
"""Q2: the fourth-order error density of a cubic, in closed form, against
measured errors of the tangent-intersection spline at equal affine-arc knots.

With p(t) = B' x B'' (a quadratic in t), affine arc length L = int |p|^(1/3) dt,
affine curvature  k = (12 p p'' - 5 p'^2) / (9 p^(8/3)),  and kappa^(1/3) =
|p|^(1/3) / |B'|, a piece of affine length h = L / n has

    E0(t) = kappa^(1/3) |k| h^4 / 8 = h^4 |12 p p'' - 5 p'^2| / (72 |p|^(7/3) |B'|)

and the one-sided tangent-intersection spline errs by max_t E0(t) / 16; the
balanced optimum by c(n) max E0 (q2_fourth_order.py), if E0 varied little.
"""
import numpy as np
import chd, schemes
from canon import SHAPES
from common import power_basis


def density(c, n, grid=2001):
    b0, b1, b2, b3 = power_basis(np.asarray(c, dtype=float))
    t = np.linspace(0, 1, grid)[:, None]
    D = b1 + t * (2 * b2 + 3 * t * b3)
    DD = 2 * b2 + 6 * t * b3
    DDD = 6 * b3
    p = D[:, 0] * DD[:, 1] - D[:, 1] * DD[:, 0]
    p1 = D[:, 0] * DDD[1] - D[:, 1] * DDD[0]                 # p'  = B' x B'''
    p2 = DD[:, 0] * DDD[1] - DD[:, 1] * DDD[0]               # p'' = B'' x B'''
    speed = np.hypot(D[:, 0], D[:, 1])
    w = np.cbrt(np.abs(p))
    L = np.sum((w[1:] + w[:-1]) / 2) / (grid - 1)
    h = L / n
    E0 = h ** 4 * np.abs(12 * p * p2 - 5 * p1 ** 2) / (72 * np.abs(p) ** (7 / 3) * speed)
    return E0


if __name__ == '__main__':
    print('T at equal affine-arc knots: measured error against max_t E0(t) / 16')
    for name in ('quarter circle R=500', 'corpus bowl arc', 'superellipse-like quarter', 'hook (line into curve)', 'flat arc, long handles'):
        c = np.asarray(SHAPES[name], dtype=float); C = c[None]
        row = []
        for n in (2, 3, 4, 6, 8):
            Q = schemes.tangents(C, schemes.affine_knots(C, n))[0]
            h = chd.hd1(c, Q, m=64)[0]
            E0 = density(c, n)
            row.append(f'n={n}: {h:.5f} / {E0.max() / 16:.5f} (density min/max {E0.min() / E0.max():.2f})')
        print(f'  {name:28s} ' + '   '.join(row))
