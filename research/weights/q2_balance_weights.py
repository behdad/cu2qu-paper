#!/usr/bin/env python3
"""Q2/Q3: universal displacements that balance a tangent-intersection spline.

Fourth-order model (see q2_fourth_order.py).  The T spline at knots of equal
affine arc length errs by  e_T(u) = u^2 (1-u)^2  on every piece, in units of
E0 = 16 x (largest one-sided error of a piece), always to the outside.  Moving
off-curve point k by  beta_k E0  toward the curve and by  a_k  forward along it
changes the error at spline parameter xi by

    - b_k(xi) [ beta_k + alpha_k (xi - k - 1/2) ],   alpha_k = - a_k theta / E0

(b_k the weight of off-curve point k in the spline, theta the turning angle of
a piece).  Minimise  max | e_T - sum_k b_k (beta_k + alpha_k (xi - k - 1/2)) |.

End tangents kept: off-curve point 0 may only slide along the start tangent,
beta_0 = alpha_0 / 2; likewise beta_{n-1} = - alpha_{n-1} / 2.

Variants:  full    beta and alpha free
           normal  alpha = 0 (move across the curve only; the two end points
                   still slide when tangents are kept)
"""
import numpy as np
from scipy.optimize import linprog

HIGHS = {'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10}


def basis(n, k, u):
    """weights of off-curve points on piece k at u: dict index -> weight"""
    w = {}
    def add(i, v):
        w[i] = w.get(i, 0.0) + v
    add(k, 2 * u * (1 - u))
    if k >= 1:
        add(k - 1, (1 - u) ** 2 / 2); add(k, (1 - u) ** 2 / 2)
    if k <= n - 2:
        add(k, u * u / 2); add(k + 1, u * u / 2)
    return w


def solve(n, tangent, variant='full', grid=129, secondary=True):
    """Returns (value, beta (n,), alpha (n,))."""
    nv = 2 * n + 1                      # beta_0..beta_{n-1}, alpha_0..alpha_{n-1}, t
    us = np.linspace(0, 1, grid)
    rows = []; rhs = []
    for k in range(n):
        for u in us:
            r = np.zeros(nv)
            for i, wgt in basis(n, k, u).items():
                r[i] += wgt
                r[n + i] += wgt * (k + u - i - 0.5)
            rows.append(r); rhs.append(u * u * (1 - u) ** 2)
    R = np.array(rows); eT = np.array(rhs)
    # |eT - R x| <= t
    A = np.vstack([np.hstack([-R[:, :2 * n], -np.ones((len(R), 1))]), np.hstack([R[:, :2 * n], -np.ones((len(R), 1))])])
    b = np.concatenate([-eT, eT])
    Aeq = []; beq = []
    if tangent:
        r = np.zeros(nv); r[0] = 1; r[n] = -0.5; Aeq.append(r); beq.append(0.0)
        r = np.zeros(nv); r[n - 1] += 1; r[2 * n - 1] += 0.5; Aeq.append(r); beq.append(0.0)
    if variant == 'normal':
        for k in range(n):
            if tangent and k in (0, n - 1):
                continue
            r = np.zeros(nv); r[n + k] = 1; Aeq.append(r); beq.append(0.0)
    cost = np.zeros(nv); cost[-1] = 1
    kw = dict(A_eq=np.array(Aeq), b_eq=np.array(beq)) if Aeq else {}
    res = linprog(cost, A_ub=A, b_ub=b, bounds=[(None, None)] * (2 * n) + [(0, None)], method='highs', options=HIGHS, **kw)
    assert res.success, res.message
    t = res.x[-1]; x = res.x[:2 * n]
    if secondary:
        # among optimal solutions: least sum |alpha| + 0.01 sum |beta - beta_mean|-ish  -> least sum |alpha|, then smooth beta
        m = 2 * n
        A2 = np.vstack([np.hstack([A[:, :2 * n], np.zeros((len(A), m))])])
        b2 = b + t * (1 + 1e-9) + 1e-12
        extra = []; eb = []
        for j in range(2 * n):
            r = np.zeros(2 * n + m); r[j] = 1; r[2 * n + j] = -1; extra.append(r); eb.append(0.0 if j >= n else 1 / 32)
            r = np.zeros(2 * n + m); r[j] = -1; r[2 * n + j] = -1; extra.append(r); eb.append(0.0 if j >= n else -1 / 32)
        A2 = np.vstack([A2, np.array(extra)]); b2 = np.concatenate([b2, eb])
        cost2 = np.concatenate([np.zeros(2 * n), 0.01 * np.ones(n), np.ones(n)])
        kw2 = dict(A_eq=np.hstack([np.array(Aeq)[:, :2 * n], np.zeros((len(Aeq), m))]), b_eq=np.array(beq)) if Aeq else {}
        res2 = linprog(cost2, A_ub=A2, b_ub=b2, bounds=[(None, None)] * (2 * n) + [(0, None)] * m, method='highs', options=HIGHS, **kw2)
        if res2.success:
            x = res2.x[:2 * n]
    val = np.abs(eT - R[:, :2 * n] @ x).max()
    return val, x[:n], x[n:]


if __name__ == '__main__':
    from q2_fourth_order import solve as solve_c1
    for tangent in (True, False):
        print('==', 'end tangents kept' if tangent else 'tangents free', '==')
        print('   n  variant   error/E0   (x16)    C1-cubic optimum   beta_k (toward the curve, units of E0)  |  alpha_k')
        for n in range(1, 9):
            ref, _ = solve_c1(n, tangent)
            for variant in ('full', 'normal'):
                v, beta, alpha = solve(n, tangent, variant)
                print(f'  {n:2d}  {variant:7s}  {v:.6f}  {16 * v:.4f}   {ref:.6f}          '
                      + ' '.join(f'{x:+.4f}' for x in beta) + '  |  ' + ' '.join(f'{x:+.4f}' for x in alpha))
        print()
