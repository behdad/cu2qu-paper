#!/usr/bin/env python3
"""Q1: the end tangents relaxed.  s_0 = x, s_{n-1} = -x with 0 <= x <= 1.

x = 1 keeps the tangents, x = 0 is the tangent-free optimum's end value.
The first off-curve point is then  c0 + (3 / 2n) (c1 - c0) - (1 - x) d3 / (4 n^3),
so the spline leaves c0 along  (c1 - c0) - (1 - x) d3 / (6 n^2).

Candidate closed form (the weights +x, -x, +x, ... of even n):
    V(x) = max_u u (1-u) (u - h),  h = (1 - x) / 2
         = (2 w^3 + 2 - 3 h - 3 h^2 + 2 h^3) / 27,  w = sqrt(1 - h + h^2).
"""
import math
import numpy as np
from q1_lp import solve
from common import phi_max


def V(x):
    h = (1 - x) / 2
    w = math.sqrt(1 - h + h * h)
    return (2 * w ** 3 + 2 - 3 * h - 3 * h * h + 2 * h ** 3) / 27


if __name__ == '__main__':
    xs = [1.0, 0.95, 0.9, 0.8, 0.75, 0.6, 0.5, 0.4, 0.25, 0.1, 0.0]
    print('least max|phi| with s_0 = x, s_{n-1} = -x  (LP), against V(x)')
    print('    x     V(x)      ' + ''.join(f'n={n:<8d}' for n in range(1, 10)))
    for x in xs:
        row = []
        for n in range(1, 10):
            if n == 1:
                # one piece: s_0 = x and s_0 = -x at once only for x = 0
                row.append('   -     ' if x != 0 else f'{phi_max([0.0]):.6f} ')
                continue
            t, s = solve(n, f'relax:{x}')
            row.append(f'{t:.6f} ')
        print(f'  {x:4.2f}  {V(x):.6f}   ' + ' '.join(row))
    print()
    print('alternating +x, -x, ... for even n has max|phi| = V(x):',
          all(abs(phi_max([x if k % 2 == 0 else -x for k in range(n)]) - V(x)) < 1e-12
              for x in xs for n in (2, 4, 6, 8)))
    # the lower bound V(|s_0|) needs only ONE pinned end weight: LP with s_0 = x and the rest free
    from scipy.optimize import linprog
    import q1_lp
    worst = 0.0
    for n in range(1, 10):
        for x in (1.0, 0.75, 0.5, 0.25, 0.0, -0.5, -1.0, 1.5):
            orig = q1_lp.solve
            # reuse solve() with a one-sided pin: emulate through the 'relax' mode on a doubled problem is
            # not needed; pin directly
            pts = [(k, u) for k in range(n) for u in np.linspace(0, 1, 65)]
            A = []; b = []
            for k, u in pts:
                c, d = q1_lp.row(n, k, u)
                A.append(np.append(c, -1.0)); b.append(-d)
                A.append(np.append(-c, -1.0)); b.append(d)
            cost = np.zeros(n + 1); cost[-1] = 1.0
            bounds = [(x, x)] + [(None, None)] * (n - 1) + [(0, None)]
            res = linprog(cost, A_ub=np.array(A), b_ub=np.array(b), bounds=bounds, method='highs', options=q1_lp.HIGHS)
            worst = max(worst, abs(phi_max(res.x[:-1]) - V(abs(x))))
    print(f'only s_0 = x pinned, the rest free, n = 1..9, x in {{1, .75, .5, .25, 0, -.5, -1, 1.5}}: |optimum - V(|x|)| <= {worst:.1e}')
    print()
    print('optimal weights for a few cases:')
    for n, x in ((3, 0.5), (4, 0.5), (5, 0.5), (6, 0.5), (4, 0.9), (5, 0.9), (6, 0.9), (7, 0.9)):
        t, s = solve(n, f'relax:{x}', antisym=True, secondary='l1')
        print(f'   n = {n}, x = {x}:  {t:.6f}  (V = {V(x):.6f})   ' + ' '.join(f'{v:+.4f}' for v in s))
    print()
    print('V(x) and the equivalent gain in n, (V(1) / V(x))^(1/3):')
    for x in xs:
        print(f'   x = {x:4.2f}   V = {V(x):.6f}   V/V(1) = {V(x) / V(1):.4f}   n shrinks by {(V(x) / V(1)) ** (1 / 3):.4f}')
    print('   slope at x = 1: dV/dx = 1/9  (V ~ 4/27 - (1 - x)/9)')

    # the angle the start tangent turns by: the spline leaves c0 along e0 - (1 - x) d3 / (6 n^2), e0 = c1 - c0
    from common import family_offcurves
    rng = np.random.default_rng(2)
    worst = 0.0
    for _ in range(2000):
        c = rng.uniform(-500, 500, size=(4, 2)); n = int(rng.integers(2, 9)); x = rng.uniform(0, 1)
        s = np.array([x if k % 2 == 0 else -x for k in range(n)])
        Q0 = family_offcurves(c, s)[0]
        e0 = c[1] - c[0]; d3 = c[3] - 3 * c[2] + 3 * c[1] - c[0]
        v = Q0 - c[0]
        ang = math.atan2(abs(e0[0] * v[1] - e0[1] * v[0]), e0 @ v)
        pred = math.atan2((1 - x) * abs(d3[0] * e0[1] - d3[1] * e0[0]) / (6 * n * n), e0 @ e0 - (1 - x) * (d3 @ e0) / (6 * n * n))
        worst = max(worst, abs(ang - pred))
    print()
    print(f'turn of the start tangent = atan2((1-x) |d3 x e0| / 6n^2, |e0|^2 - (1-x) d3.e0 / 6n^2): worst gap on 2000 random cubics {worst:.1e} rad')
