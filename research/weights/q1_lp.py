#!/usr/bin/env python3
"""Q1.  Best weights for uniform splits under the equal-parameter error.

Independent of the earlier weights_lp.py: the linear program is solved by an
exchange (cutting-plane) loop whose constraints are the exact extrema of each
piece's cubic phi_k, not a fixed grid, so the optimum is exact to ~1e-12.

    minimise  t   subject to  |phi_k(u; s)| <= t  for all pieces k and u in [0, 1]

phi_k(u) = (1-u)^2 g_k + u^2 g_{k+1} + u (1-u) (u - h_k) is linear in s.

Modes (which weights are pinned):
  both      s_0 = +1, s_{n-1} = -1        both end tangents kept
  one       s_0 = +1                      only the start tangent kept
  free      nothing                       no tangent kept
  relax:x   s_0 = x, s_{n-1} = -x         end tangents off by (1-x) d / 4
"""
import sys
import numpy as np
from scipy.optimize import linprog
from common import phi_coeffs, cubic_extrema, poly, phi_max, ramp, alternating


HIGHS = {'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10}


def row(n, k, u):
    """phi_k(u) = c . s + const."""
    c = np.zeros(n)
    if k >= 1:
        c[k - 1] += (1 - u) ** 2 / 8
        c[k] += (1 - u) ** 2 / 8
    if k <= n - 2:
        c[k] += u * u / 8
        c[k + 1] += u * u / 8
    c[k] += u * (1 - u) / 2
    return c, u * (1 - u) * (u - 0.5)


def solve(n, mode='both', antisym=False, tol=1e-12, verbose=False, secondary=None):
    """Returns (t, s).  antisym: additionally impose s_k = -s_{n-1-k}.
    secondary: None, or 'l1' to pick, among optimal s, one of least sum |s_k|
    (makes the reported weights canonical when the optimum is not unique)."""
    pts = [(k, u) for k in range(n) for u in np.linspace(0, 1, 17)]
    bounds = [(None, None)] * n + [(0, None)]
    if mode == 'both':
        bounds[0] = (1.0, 1.0); bounds[n - 1] = (-1.0, -1.0)
    elif mode == 'one':
        bounds[0] = (1.0, 1.0)
    elif mode.startswith('relax:'):
        x = float(mode.split(':')[1])
        bounds[0] = (x, x); bounds[n - 1] = (-x, -x)
    elif mode != 'free':
        raise ValueError(mode)
    Aeq = []; beq = []
    if antisym:
        for k in range((n + 1) // 2):
            r = np.zeros(n + 1); r[k] += 1; r[n - 1 - k] += 1
            Aeq.append(r); beq.append(0.0)
    for it in range(200):
        A = []; b = []
        for k, u in pts:
            c, d = row(n, k, u)
            A.append(np.append(c, -1.0)); b.append(-d)
            A.append(np.append(-c, -1.0)); b.append(d)
        cost = np.zeros(n + 1); cost[-1] = 1.0
        res = linprog(cost, A_ub=np.array(A), b_ub=np.array(b),
                      A_eq=np.array(Aeq) if Aeq else None, b_eq=np.array(beq) if Aeq else None,
                      bounds=bounds, method='highs', options=HIGHS)
        assert res.success, res.message
        t, s = res.x[-1], res.x[:-1]
        co = phi_coeffs(s)
        worst = 0.0; new = []
        for k in range(n):
            for u in cubic_extrema(co[k]):
                v = abs(poly(co[k], u))
                worst = max(worst, v)
                if v > t + tol:
                    new.append((k, u))
        if verbose:
            print(f'   it {it}: lp {t:.15f} true {worst:.15f} added {len(new)}')
        if not new:
            break
        pts.extend(new)
    if secondary == 'l1':
        # among s with max |phi| <= t (on the collected points), least sum |s_k|
        m = n
        A = []; b = []
        for k, u in pts:
            c, d = row(n, k, u)
            A.append(np.concatenate([c, np.zeros(m)])); b.append(t + 1e-12 - d)
            A.append(np.concatenate([-c, np.zeros(m)])); b.append(t + 1e-12 + d)
        for k in range(n):
            r = np.zeros(n + m); r[k] = 1; r[n + k] = -1; A.append(r); b.append(0.0)
            r = np.zeros(n + m); r[k] = -1; r[n + k] = -1; A.append(r); b.append(0.0)
        cost = np.concatenate([np.zeros(n), np.ones(m)])
        Aeq2 = [np.concatenate([r[:n], np.zeros(m)]) for r in Aeq]
        res = linprog(cost, A_ub=np.array(A), b_ub=np.array(b),
                      A_eq=np.array(Aeq2) if Aeq2 else None, b_eq=np.array(beq) if Aeq2 else None,
                      bounds=bounds[:n] + [(0, None)] * m, method='highs', options=HIGHS)
        if res.success:
            s2 = res.x[:n]
            if phi_max(s2) <= t + 1e-9:
                s = s2
    return phi_max(s), s


def fmt(s, digits=4, limit=11):
    s = list(s)
    txt = ' '.join(f'{x:+.{digits}f}' for x in s[:limit])
    return txt + (' ...' if len(s) > limit else '')


if __name__ == '__main__':
    print('== both end tangents kept (s_0 = +1, s_{n-1} = -1) ==')
    print('   n   optimum Phi        ramp Phi    alternating   optimum/(4/27)   n^3/Phi(opt)  n^3/Phi(ramp)   weights (antisymmetric, least sum |s|)')
    for n in list(range(2, 17)) + [20, 21, 30, 31]:
        t, s = solve(n, 'both', antisym=True, secondary='l1')
        t2, _ = solve(n, 'both', antisym=False)
        assert abs(t - t2) < 2e-9, (n, t, t2)       # symmetrising loses nothing
        alt = f'{phi_max(alternating(n)):.6f}' if n % 2 == 0 else '    -   '
        print(f'{n:4d}   {t:.12f}   {phi_max(ramp(n)):.6f}    {alt}      {t / (4 / 27):.6f}    '
              f'{n ** 3 / t:12.2f}  {n ** 3 / phi_max(ramp(n)):12.2f}    {fmt(s)}')
    print('   4/27 =', 4 / 27, '  27/128 =', 27 / 128)

    print()
    print('== only the start tangent kept (s_0 = +1, the rest free) ==')
    for n in range(1, 13):
        t, s = solve(n, 'one')
        print(f'{n:4d}   {t:.12f}   {fmt(s)}')

    print()
    print('== no tangent kept ==')
    for n in (1, 2, 3, 4, 5, 8, 13):
        t, s = solve(n, 'free')
        print(f'{n:4d}   {t:.12f}   sqrt(3)/36 = {3 ** 0.5 / 36:.12f}   {fmt(s)}')
