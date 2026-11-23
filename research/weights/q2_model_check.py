#!/usr/bin/env python3
"""Check the fourth-order model (q2_fourth_order.py) against numerically
optimised splines on cubic arcs of a circle.

For a circle of radius R and a piece turning by theta, E0 = R theta^4 / 8, so
the model predicts, for an arc of total angle Theta cut into n:

    best distance = c(n) R (Theta / n)^4 / 8,

c(n) the model constant (tangents kept or free).  The one-sided tangent-
intersection quadratic has c = 1/16, exactly R (1 - cos(theta/2))^2 / (2 cos(theta/2)).
"""
import math
import numpy as np
import chd, schemes
from opt import best_spline
from q2_fourth_order import solve


def arc_cubic(R, Theta):
    k = 4 / 3 * math.tan(Theta / 4)
    c0 = np.array([R, 0.0]); c3 = R * np.array([math.cos(Theta), math.sin(Theta)])
    c1 = c0 + k * R * np.array([0.0, 1.0])
    c2 = c3 - k * R * np.array([-math.sin(Theta), math.cos(Theta)])
    return np.array([c0, c1, c2, c3])


if __name__ == '__main__':
    R = 1000.0
    for deg in (20, 45, 90):
        Theta = math.radians(deg)
        c = arc_cubic(R, Theta); C = c[None]
        print(f'circular arc of {deg} degrees, R = {R:g}  (the cubic itself is within {R * 4 / 27 * math.sin(Theta / 4) ** 6 / math.cos(Theta / 4) ** 2:.2e} of the circle)')
        print('   n    T (measured)  /E0     model 1/16 |  TAN* (measured)  /E0    model   |  FREE* (measured)  /E0    model')
        for n in range(1, 7):
            E0 = R * (Theta / n) ** 4 / 8
            Tu = schemes.tangents(C, schemes.uniform_knots(1, n))[0]
            Ta = schemes.tangents(C, schemes.affine_knots(C, n))[0]
            hT = chd.hd1(c, Ta, m=64)[0]
            starts = [Ta, Tu] + ([schemes.optu(C, n)[0]] if n > 1 else [])
            ht = best_spline(c, n, 'tan', starts)[0] if n > 1 else hT
            hf = best_spline(c, n, 'free', starts + [schemes.zero(C, n)[0]])[0]
            mt, _ = solve(n, True); mf, _ = solve(n, False)
            print(f'  {n:2d}   {hT:12.6f}  {hT / E0:.5f}  {1 / 16:.5f}  |  {ht:12.6f}  {ht / E0:.5f}  {mt:.5f}  |  {hf:12.6f}  {hf / E0:.5f}  {mf:.5f}')
        print()
