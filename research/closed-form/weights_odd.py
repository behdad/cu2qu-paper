#!/usr/bin/env python3
"""Optimal tangent-keeping weights for odd n, and how each n's threshold
(the largest |d3| / tolerance it can take) compares with cu2qu's ramp."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from weights_lp import solve, phi_max

print('  n   optimal  ramp     threshold n^3/Phi: optimal   ramp     optimal weights')
prev = 0
for n in range(2, 16):
    ramp = np.array([1 - 2 * k / (n - 1) for k in range(n)])
    t, s = solve(n, 'tangent')
    w = ' '.join(f'{x:+.2f}' for x in s)
    print(f'{n:4d}  {t:.4f}   {phi_max(ramp):.4f}   {n**3 / t:18.1f}  {n**3 / phi_max(ramp):9.1f}     {w}')
for n in (17, 21, 31):
    t, s = solve(n, 'tangent')
    print(f'{n:4d}  {t:.4f}   first weights: ' + ' '.join(f'{x:+.2f}' for x in s[:8]) + ' ... middle ' + f'{s[n // 2]:+.2f}')
