#!/usr/bin/env python3
"""Q2/Q3: third-order constants of the geometric error for EQUAL-parameter splits.

For a spline built on n equal pieces whose on-curve points stay within O(l^2)
of the split points, the geometric error is, to leading order,

    c * |d3 x T| / n^3      (T the unit tangent where the error peaks),

only the part of d3 across the curve counts.  Predicted c:
    cu2qu ramp            Phi(n)   (4/27, 27/128, 0.2353, ...)
    alternating           4/27
    T at equal parameters 2/27 in the two end pieces, sqrt(3)/36 inside
    s = 0                 sqrt(3)/36

Test curve: y = a x^2 + b x^3 over a long base (d3 = (0, b) is across the
curve; b << a so that the fourth-order term, which is of order b^2/a, is
negligible), measured Hausdorff distance times n^3 / b.
"""
import numpy as np
import chd, schemes
from common import phi_max, ramp as ramp_w

L, a, b = 1000.0, 40.0, 4.0
c = np.array([[0, 0], [L / 3, 0], [2 * L / 3, a / 3], [L, a + b]], dtype=float)
C = c[None]
print(f'cubic (L t, {a:g} t^2 + {b:g} t^3), L = {L:g}; measured Hausdorff distance * n^3 / |d3|')
print('   n     ramp    Phi(n)  |  alternating/optu  Phi_opt  |  T equal-t   2/27    |  s = 0    sqrt(3)/36')
for n in range(2, 9):
    row = []
    for Q in (schemes.ramp(C, n)[0], schemes.optu(C, n)[0], schemes.tangents(C, schemes.uniform_knots(1, n))[0], schemes.zero(C, n)[0]):
        row.append(chd.hd1(c, Q, m=64)[0] * n ** 3 / b)
    print(f'  {n:2d}   {row[0]:.4f}  {phi_max(ramp_w(n)):.4f}  |      {row[1]:.4f}       {phi_max(schemes.optu_weights(n)):.4f}   |  {row[2]:.4f}    {2 / 27:.4f}  |  {row[3]:.4f}   {3 ** 0.5 / 36:.4f}')
