#!/usr/bin/env python3
"""Two quick facts about splines of n quadratics over equal pieces of a cubic.

Write the off-curve point of piece k as  q_M + s_k * e / 4,  where q_M is
the control point of least equal-parameter error, (3 (b + c) - a - d) / 4,
e = d3 / n^3, and on-curve points are implied at the midpoints.
  cu2qu:               s_k = 1 - 2k/(n-1)   (a linear ramp from +1 to -1)
  closed-form split:   s_k = 0              (kurbo, lyon)

1. With s_k = 0, are the split points of the cubic exactly the midpoints of
   adjacent off-curve points?  (Then that output, too, has implied on-curve
   points in TrueType.)
2. What is the largest error, as a multiple of |d3| / n^3, for a few other
   choices of s that keep the end tangents (s_0 = 1, s_{n-1} = -1)?
"""
import math, random

def pieces(c, n):
    c0, c1, c2, c3 = c
    b1 = 3 * (c1 - c0); b2 = 3 * (c0 - 2 * c1 + c2); b3 = c3 - 3 * c2 + 3 * c1 - c0
    B = lambda t: c0 + t * (b1 + t * (b2 + t * b3))
    D = lambda t: b1 + t * (2 * b2 + 3 * t * b3)
    out = []
    for k in range(n):
        t0, t1 = k / n, (k + 1) / n
        a, d = B(t0), B(t1)
        out.append((a, a + D(t0) / (3 * n), d - D(t1) / (3 * n), d))
    return out, b3

random.seed(3)
worst = 0.0
for _ in range(2000):
    c = [complex(random.uniform(-1000, 1000), random.uniform(-1000, 1000)) for _ in range(4)]
    n = random.randint(2, 30)
    ps, d3 = pieces(c, n)
    off = [(3 * (b + cc) - a - d) / 4 for a, b, cc, d in ps]
    for k in range(1, n):
        mid = (off[k - 1] + off[k]) / 2
        worst = max(worst, abs(mid - ps[k][0]) / max(1.0, abs(d3)))
print(f'1. closed-form split: split point vs midpoint of adjacent off-curve points, largest relative gap {worst:.2e}')

def Phi_for(s):
    """Largest |phi_k(u)| for weights s (list of n), by dense evaluation."""
    n = len(s)
    g = [0.0] + [(s[k - 1] + s[k]) / 8.0 for k in range(1, n)] + [0.0]
    best = 0.0
    for k in range(n):
        h = (1 - s[k]) / 2.0
        co = (g[k], -(2 * g[k] + h), 1 + g[k] + g[k + 1] + h, -1.0)
        for j in range(0, 2001):
            u = j / 2000.0
            best = max(best, abs(co[0] + u * (co[1] + u * (co[2] + u * co[3]))))
    return best

print('2. largest error / (|d3| / n^3):')
print('     n   cu2qu ramp   +1,0..0,-1   alternating   s=0 (tangents not kept)')
for n in (2, 3, 4, 6, 8, 12, 20, 40):
    ramp = [1 - 2 * k / (n - 1) for k in range(n)]
    ends = [1.0] + [0.0] * (n - 2) + [-1.0] if n > 2 else [1.0, -1.0]
    alt = [1.0 if k % 2 == 0 else -1.0 for k in range(n)]
    alt_txt = f'{Phi_for(alt):.4f}' if n % 2 == 0 else '  (odd n)'
    print(f'  {n:4d}    {Phi_for(ramp):.4f}       {Phi_for(ends):.4f}       {alt_txt}        {Phi_for([0.0] * n):.4f}')
