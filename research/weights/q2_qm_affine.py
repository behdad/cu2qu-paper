#!/usr/bin/env python3
"""Q2: the "catch" of unequal pieces with implied midpoints, measured.

Take the pieces between knots of equal affine arc length and give each its own
d3-family point with s = 0, q_M = (B(t0) + B(t1))/2 - l (B'(t1) - B'(t0))/4.
The implied on-curve points are then off the split points by about
(l_k - l_{k-1})/4 B'(t_k) -- along the curve -- so the distance at equal
(piecewise) parameters is large, of first order.  The distance between the
curves is nevertheless of fourth order: the same order as the
tangent-intersection spline at those knots, and far below s = 0 on equal pieces.
"""
import numpy as np
import chd, schemes
from canon import SHAPES
from common import power_basis


def qm_at(C, knots):
    b0, b1, b2, b3 = power_basis(C)
    t = knots[:, :, None]
    P = b0[:, None] + t * (b1[:, None] + t * (b2[:, None] + t * b3[:, None]))
    D = b1[:, None] + t * (2 * b2[:, None] + 3 * t * b3[:, None])
    l = (knots[:, 1:] - knots[:, :-1])[..., None]
    return (P[:, :-1] + P[:, 1:]) / 2 - l * (D[:, 1:] - D[:, :-1]) / 4


def piecewise_equal_param(c, Q, knots, m=64):
    """largest |S(k + u) - B(t_k + u l_k)|"""
    b0, b1, b2, b3 = power_basis(c)
    n = len(Q)
    on = [c[0]] + [(Q[k - 1] + Q[k]) / 2 for k in range(1, n)] + [c[3]]
    worst = 0.0
    for k in range(n):
        u = np.linspace(0, 1, m)[:, None]
        t = knots[k] + u * (knots[k + 1] - knots[k])
        B = b0 + t * (b1 + t * (b2 + t * b3))
        S = (1 - u) ** 2 * on[k] + 2 * u * (1 - u) * Q[k] + u * u * on[k + 1]
        worst = max(worst, np.hypot(*(S - B).T).max())
    return worst


if __name__ == '__main__':
    for name in ('corpus bowl arc', 'hook (line into curve)', 'flat, one-sided handles'):
        c = np.asarray(SHAPES[name], dtype=float); C = c[None]
        print(name)
        print('    n   q_M at affine knots: Hausdorff   at equal piecewise parameter |  T at affine knots  |  s = 0 on equal pieces')
        for n in (2, 3, 4, 6, 8):
            ka = schemes.affine_knots(C, n)
            Qa = qm_at(C, ka)[0]
            print(f'   {n:2d}        {chd.hd1(c, Qa, m=64)[0]:10.5f}                  {piecewise_equal_param(c, Qa, ka[0]):10.4f}'
                  f'           |    {chd.hd1(c, schemes.tangents(C, ka)[0], m=64)[0]:10.5f}       |    {chd.hd1(c, schemes.zero(C, n)[0], m=64)[0]:10.5f}')
