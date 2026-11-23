#!/usr/bin/env python3
"""Best weights for the uniform quadratic spline of a cubic, per n.

Family.  Cut the cubic into n equal pieces; piece k gets the off-curve
point  q_M(k) + s_k * e / 4,  with q_M the control point of least
equal-parameter error for the piece, e = d3 / n^3, and on-curve points
implied halfway between consecutive off-curve points.  cu2qu is the ramp
s_k = 1 - 2k/(n-1); the closed-form split is s = 0.

The difference between quadratic k and its piece at parameter u is
e * phi_k(u; s) with

    phi_k = g_k - u (2 g_k + h_k) + u^2 (1 + g_k + g_{k+1} + h_k) - u^3
    h_k = (1 - s_k) / 2,   g_k = (s_{k-1} + s_k) / 8,   g_0 = g_n = 0

which is linear in s.  So minimizing the largest |phi| is a linear program.

Constraint sets:
  'tangent': s_0 = +1 and s_{n-1} = -1 (the first and last off-curve
             points on the cubic's end tangents), the rest free.
  'free':    all free (tangents not kept); expect s = 0, 0.0481.
"""
import math, sys
import numpy as np
from scipy.optimize import linprog

U = np.linspace(0.0, 1.0, 401)

def rows(n):
    """For each (k, u): coefficients c (n,) and constant d with phi = c.s + d."""
    R = []
    for k in range(n):
        for u in U:
            c = np.zeros(n)
            d = 0.0
            # g_k = (s_{k-1} + s_k)/8 for 1 <= k <= n-1
            def add_g(idx, coef):
                if 1 <= idx <= n - 1:
                    c[idx - 1] += coef / 8.0
                    c[idx] += coef / 8.0
            # phi = g_k (1 - 2u + u^2) + g_{k+1} u^2 + h_k (u^2 - u) + u^2 - u^3
            add_g(k, (1 - u) ** 2)
            add_g(k + 1, u * u)
            # h_k = (1 - s_k)/2
            d += 0.5 * (u * u - u)
            c[k] += -0.5 * (u * u - u)
            d += u * u - u ** 3
            R.append((c, d))
    return R

def solve(n, mode):
    R = rows(n)
    # variables: s_0..s_{n-1}, t
    A = []; b = []
    for c, d in R:
        A.append(np.append(c, -1.0)); b.append(-d)      #  c.s + d <= t
        A.append(np.append(-c, -1.0)); b.append(d)      # -c.s - d <= t
    cost = np.zeros(n + 1); cost[-1] = 1.0
    bounds = [(None, None)] * n + [(0, None)]
    if mode == 'tangent':
        bounds[0] = (1.0, 1.0)
        bounds[n - 1] = (-1.0, -1.0)
    res = linprog(cost, A_ub=np.array(A), b_ub=np.array(b), bounds=bounds, method='highs')
    assert res.success, res.message
    return res.x[-1], res.x[:-1]

def phi_max(s):
    n = len(s)
    best = 0.0
    for c, d in rows(n):
        best = max(best, abs(float(np.dot(c, s)) + d))
    return best

if __name__ == '__main__':
    ns = [2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 24, 32]
    print('  n   ramp     ends-only  alternating  optimal(tangent)  gain in n vs ramp   optimal weights')
    for n in ns:
        ramp = np.array([1 - 2 * k / (n - 1) for k in range(n)])
        ends = np.array([1.0] + [0.0] * (n - 2) + [-1.0]) if n > 2 else np.array([1.0, -1.0])
        alt = np.array([1.0 if k % 2 == 0 else -1.0 for k in range(n)])
        t, s = solve(n, 'tangent')
        alt_txt = f'{phi_max(alt):.4f}' if n % 2 == 0 else '   -  '
        w = ' '.join(f'{x:+.3f}' for x in s[:min(n, 8)]) + (' ...' if n > 8 else '')
        print(f'{n:4d}  {phi_max(ramp):.4f}   {phi_max(ends):.4f}     {alt_txt}       {t:.4f}          '
              f'{(phi_max(ramp) / t) ** (1 / 3):.3f}x            {w}')
    t, s = solve(6, 'free')
    print('free (tangents not kept), n=6:', f'{t:.4f}', ' '.join(f'{x:+.3f}' for x in s))
