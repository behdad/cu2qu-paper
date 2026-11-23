#!/usr/bin/env python3
"""Q3 under the equal-parameter measure: let the two end off-curve points slide
along the end tangents (and every interior point move freely in the plane).

Uniform split.  Write Q_k = q_M,k + (|d|/4) (s_k dh + r_k dp), dh the unit
vector of d3 and dp the one across it.  The error at equal parameters is then
|d| (phi_k(u; s) dh + omega_k(u; r) dp), omega the linear part of phi.  The
start tangent is kept iff (s_0 - 1) dh + r_0 dp is along the tangent, i.e.

    r_0 = (s_0 - 1) tan psi_0,      r_{n-1} = (s_{n-1} + 1) tan psi_1,

psi the angle from d3 to the end tangent.  The d3 family is r = 0, s_0 = 1,
s_{n-1} = -1.  Minimise  max sqrt(phi^2 + omega^2)  (a second-order cone
program; solved as a linear program with the disc replaced by a 128-gon,
accurate to 3e-4 relative).

So the constant depends on the cubic only through psi_0 and psi_1.
"""
import math, sys
import numpy as np
from scipy.optimize import linprog
from q1_lp import row
from common import CORPUS

HIGHS = {'primal_feasibility_tolerance': 1e-9, 'dual_feasibility_tolerance': 1e-9}


def solve(n, tan0, tan1, grid=25, ndir=128):
    """min max |error| / |d| with the end points on the end tangents.  Returns (value, s, r)."""
    # variables: s (n), r (n), t
    nv = 2 * n + 1
    th = np.linspace(0, 2 * math.pi, ndir, endpoint=False)
    ct = np.cos(th); st = np.sin(th)
    A = []; b = []
    for k in range(n):
        for u in np.linspace(0, 1, grid):
            c, const = row(n, k, u)
            # phi = c.s + const ; omega = c.r ;  phi ct + omega st <= t
            blk = np.zeros((ndir, nv))
            blk[:, :n] = ct[:, None] * c[None, :]
            blk[:, n:2 * n] = st[:, None] * c[None, :]
            blk[:, -1] = -1.0
            A.append(blk); b.append(-const * ct)
    A = np.vstack(A); b = np.concatenate(b)
    Aeq = []; beq = []
    if math.isinf(tan0):
        r = np.zeros(nv); r[0] = 1; Aeq.append(r); beq.append(1.0)          # s_0 = 1 (d3 across the tangent: no sliding in d3)
        # r_0 free: sliding is purely across d3
    else:
        r = np.zeros(nv); r[n] = 1; r[0] = -tan0; Aeq.append(r); beq.append(-tan0)       # r_0 = (s_0 - 1) tan0
    if math.isinf(tan1):
        r = np.zeros(nv); r[n - 1] = 1; Aeq.append(r); beq.append(-1.0)
    else:
        r = np.zeros(nv); r[2 * n - 1] += 1; r[n - 1] += -tan1; Aeq.append(r); beq.append(tan1)   # r_{n-1} = (s_{n-1} + 1) tan1
    cost = np.zeros(nv); cost[-1] = 1
    res = linprog(cost, A_ub=A, b_ub=b, A_eq=np.array(Aeq), b_eq=np.array(beq),
                  bounds=[(None, None)] * (2 * n) + [(0, None)], method='highs', options=HIGHS)
    assert res.success, res.message
    s = res.x[:n]; r = res.x[n:2 * n]
    # true value on a fine grid
    worst = 0.0
    for k in range(n):
        for u in np.linspace(0, 1, 201):
            c, const = row(n, k, u)
            worst = max(worst, math.hypot(c @ s + const, c @ r))
    return worst, s, r


if __name__ == '__main__':
    print('least equal-parameter error / (|d3| / n^3), end tangents kept, end points free to slide')
    print('symmetric arc: psi_0 = psi, psi_1 = -psi (d3 along the chord)')
    degs = [0, 10, 20, 30, 45, 60, 75, 85, 90]
    print('   n   ' + ''.join(f'psi={d:<5d}' for d in degs) + '   d3 family (no sliding)')
    from q1_lp import solve as solve_d3
    for n in (1, 2, 3, 4, 5, 6):
        vals = []
        for dg in degs:
            tn = math.inf if dg == 90 else math.tan(math.radians(dg))
            if n == 1 and dg == 90:
                vals.append(float('nan')); continue
            v, s, r = solve(n, tn, -tn if not math.isinf(tn) else math.inf)
            vals.append(v)
        ref = solve_d3(n, 'both')[0] if n > 1 else float('nan')
        print(f'  {n:2d}   ' + ''.join(f'{v:<9.4f}' for v in vals) + f'   {ref:.4f}')
    print('   (4/27 = 0.1481, 27/128 = 0.2109, sqrt(3)/36 = 0.0481)')
    print()
    print('weights at psi = 45 degrees:')
    for n in (2, 3, 4):
        v, s, r = solve(n, 1.0, -1.0)
        print(f'   n = {n}: {v:.4f}   s = ' + ' '.join(f'{x:+.3f}' for x in s) + '   r = ' + ' '.join(f'{x:+.3f}' for x in r))

    # corpus: distribution of psi at the two ends
    import json
    C = np.array(json.load(open(CORPUS)), dtype=float)
    d3 = C[:, 3] - 3 * C[:, 2] + 3 * C[:, 1] - C[:, 0]
    from schemes import unit_end_tangents, cross
    e0, e1 = unit_end_tangents(C)
    L = np.hypot(d3[:, 0], d3[:, 1])
    ok = L > 1e-9
    psi0 = np.degrees(np.arcsin(np.clip(np.abs(cross(d3[ok], e0[ok])) / L[ok], 0, 1)))
    psi1 = np.degrees(np.arcsin(np.clip(np.abs(cross(d3[ok], e1[ok])) / L[ok], 0, 1)))
    print()
    print('corpus: angle between d3 and the end tangents (degrees, folded to 0..90)')
    for nm, a in (('start', psi0), ('end', psi1), ('larger of the two', np.maximum(psi0, psi1))):
        print(f'   {nm:18s} quantiles 10/25/50/75/90%: ' + ' '.join(f'{v:5.1f}' for v in np.quantile(a, [.1, .25, .5, .75, .9])))
