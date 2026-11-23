#!/usr/bin/env python3
"""Q2/Q3: the fourth-order model of the geometric error, and its universal constants.

Where the cubic has no inflection, parametrise it by affine arc length sigma
(d sigma = |B' x B''|^(1/3) dt).  Locally it is  y = x^2/2 + k x^4/8  in affine
coordinates (k the affine curvature), so against any parabola the third-order
term is gone.  An implied-on-curve quadratic spline is a C1 piecewise quadratic
in its own uniform parameter; moving its off-curve points across the curve
changes the distance by a C1 piecewise quadratic, moving them along the curve
by xi times one.  Together: any C1 piecewise cubic.  So, to leading order, with
h the affine length of a piece and E0 = k h^4 / 8 (in affine units; times
kappa^(1/3) for a Euclidean distance):

    distance / E0  =  min over C1 piecewise cubics c of  max | xi^4 - c(xi) |

on n unit pieces, with the end values zero (the spline ends on the curve) and,
when the end tangents are kept, the end slopes zero too.

Reference values: one quadratic through both end tangents (u^2 (1-u)^2) has
1/16; far from the ends the best is u^2 (1-u)^2 - 1/32, i.e. 1/32.
"""
import numpy as np
from scipy.optimize import linprog

HIGHS = {'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10}


def solve(n, tangent, grid=65, rounds=6):
    """min sup |e|, e_j(u) = u^4 + a u^3 + b u^2 + c u + d on piece j.  Returns (value, coeffs (n, 4))."""
    nv = 4 * n + 1
    def val_row(j, u):       # value of piece j at u, without the u^4 term
        r = np.zeros(nv); r[4 * j:4 * j + 4] = [u ** 3, u ** 2, u, 1.0]; return r
    def der_row(j, u):
        r = np.zeros(nv); r[4 * j:4 * j + 4] = [3 * u * u, 2 * u, 1.0, 0.0]; return r
    Aeq = []; beq = []
    for j in range(n - 1):
        Aeq.append(val_row(j, 1.0) - val_row(j + 1, 0.0)); beq.append(-1.0)          # e_j(1) = e_{j+1}(0)
        Aeq.append(der_row(j, 1.0) - der_row(j + 1, 0.0)); beq.append(-4.0)
    Aeq.append(val_row(0, 0.0)); beq.append(0.0)
    Aeq.append(val_row(n - 1, 1.0)); beq.append(-1.0)
    if tangent:
        Aeq.append(der_row(0, 0.0)); beq.append(0.0)
        Aeq.append(der_row(n - 1, 1.0)); beq.append(-4.0)
    pts = [(j, u) for j in range(n) for u in np.linspace(0, 1, grid)]
    cost = np.zeros(nv); cost[-1] = 1.0
    for rd in range(rounds):
        A = []; b = []
        for j, u in pts:
            r = val_row(j, u)
            rp = r.copy(); rp[-1] = -1.0; A.append(rp); b.append(-u ** 4)
            rm = -r; rm[-1] = -1.0; A.append(rm); b.append(u ** 4)
        res = linprog(cost, A_ub=np.array(A), b_ub=np.array(b), A_eq=np.array(Aeq), b_eq=np.array(beq),
                      bounds=[(None, None)] * (nv - 1) + [(0, None)], method='highs', options=HIGHS)
        assert res.success, res.message
        co = res.x[:-1].reshape(n, 4)
        # exact extrema of each quartic
        worst = 0.0; new = []
        for j in range(n):
            a, bb, c, d = co[j]
            roots = np.roots([4.0, 3 * a, 2 * bb, c])
            for u in [0.0, 1.0] + [r.real for r in roots if abs(r.imag) < 1e-12 and 0 < r.real < 1]:
                v = abs(u ** 4 + a * u ** 3 + bb * u * u + c * u + d)
                worst = max(worst, v)
                if v > res.x[-1] + 1e-12:
                    new.append((j, u))
        if not new:
            break
        pts.extend(new)
    return worst, co


if __name__ == '__main__':
    print('fourth-order constants: best distance / E0, E0 = k h^4 / 8, n pieces of equal affine length h')
    print('   n   tangents kept   (x16)     tangents free   (x16)')
    for n in list(range(1, 11)) + [16, 24]:
        t, _ = solve(n, True)
        f, _ = solve(n, False)
        print(f'  {n:2d}   {t:.8f}   {16 * t:.5f}     {f:.8f}   {16 * f:.5f}')
    print('   references: 1/16 =', 1 / 16, '  1/32 =', 1 / 32, '  (3 - 2 sqrt 2)/16 =', (3 - 2 * 2 ** 0.5) / 16)
    print()
    print('as a function of the whole curve (affine length L = n h): distance = C(n) k L^4 / 8 with C(n) = const / n^4')
    for n in (1, 2, 3, 4, 5, 6, 8):
        t, _ = solve(n, True); f, _ = solve(n, False)
        print(f'  n = {n}: kept {t / n ** 4:.3e}   free {f / n ** 4:.3e}   kept/free {t / f:.2f}')
