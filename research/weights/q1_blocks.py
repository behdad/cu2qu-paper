#!/usr/bin/env python3
"""Q1: the alternating spline as blocks, and simple weights for odd n.

1. Alternating weights on n = 2m equal pieces are exactly: cut the cubic into m
   equal sub-cubics and give each cu2qu's own two-quadratic spline.  Checked
   against fontTools (split_cubic_into_n_iter + cubic_approx_spline(.., 2, ..)).

2. Blocks.  Where s_{k-1} = -1 and s_k = +1 the on-curve point between them is
   on the cubic (g_k = 0) and neither side feels the other, so a weight vector
   made of blocks that each start with +1 and end with -1 has the largest of
   its blocks' constants.  Odd n >= 5 can therefore be: pairs (+1, -1) and one
   block of 5.
"""
import sys
import numpy as np
from common import phi_max, alternating, family_offcurves, FT_LIB
from schemes import optu_weights

sys.path.insert(0, FT_LIB)
from fontTools.cu2qu.cu2qu import cubic_approx_spline, split_cubic_into_n_iter

rng = np.random.default_rng(4)
worst = 0.0
for trial in range(5000):
    c = rng.uniform(-1000, 1000, size=(4, 2))
    m = int(rng.integers(1, 9))
    cc = [complex(*p) for p in c]
    pts = []
    for sub in split_cubic_into_n_iter(cc[0], cc[1], cc[2], cc[3], m):
        sp = cubic_approx_spline(sub, 2, 1e30, True)
        pts.extend(sp[1:-1])
    ft = np.array([[z.real, z.imag] for z in pts])
    Q = family_offcurves(c, alternating(2 * m))
    worst = max(worst, np.abs(ft - Q).max())
print(f'1. alternating on 2m pieces = m sub-cubics x cu2qu n=2: worst gap over 5000 random cubics, m = 1..8: {worst:.2e}')

print('2. odd n: constants of simple patterns')
a5 = optu_weights(5)
print(f'   optimal block of 5: {" ".join(f"{x:+.6f}" for x in a5)}   Phi = {phi_max(a5):.9f}')
b5 = np.array([1.0, -1.0, 0.0, 1.0, -1.0])
print(f'   (+1, -1, 0, +1, -1): Phi = {phi_max(b5):.9f}   closed form u^3 - 15u^2/8 + u at u = (15 - sqrt 33)/24: '
      f'{(lambda u: u ** 3 - 15 * u * u / 8 + u)((15 - 33 ** 0.5) / 24):.9f}')
print('     n   LP optimum   pairs + optimal 5-block   pairs + (+1,-1,0,+1,-1)   pairs + (+1,0,-1)   ramp')
from common import ramp
for n in (5, 7, 9, 11, 13, 15):
    m = (n - 5) // 2
    left = m // 2; right = m - left
    def with_block(blk):
        return np.concatenate([np.tile([1.0, -1.0], left), blk, np.tile([1.0, -1.0], right)])
    t3 = np.concatenate([np.tile([1.0, -1.0], (n - 3) // 2), [1.0, 0.0, -1.0]])
    print(f'    {n:2d}   {phi_max(optu_weights(n)):.6f}     {phi_max(with_block(a5)):.6f}                 {phi_max(with_block(b5)):.6f}'
          f'                  {phi_max(t3):.6f}            {phi_max(ramp(n)):.6f}')
print('   block position does not matter:',
      all(abs(phi_max(np.concatenate([np.tile([1.0, -1.0], i), a5, np.tile([1.0, -1.0], 4 - i)])) - phi_max(a5)) < 1e-12 for i in range(5)))
