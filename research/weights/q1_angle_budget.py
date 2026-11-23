#!/usr/bin/env python3
"""Q1/Q5: what a bounded end-tangent error buys, on the corpus, with n from the
closed-form equal-parameter bound.

Scheme with budget delta (degrees).  n >= 2: weights +x, -x, +x, ... (even n;
the LP optimum with s_0 = x, s_{n-1} = -x for odd n), x the smallest value in
[0, 1] for which neither end tangent turns by more than delta:

    the spline leaves c0 along  e0 - (1 - x) d3 / (6 n^2),  e0 = c1 - c0
    =>  1 - x <= 6 n^2 tan(delta) |e0|^2 / (|d3 x e0| + tan(delta) d3.e0)

and fits iff Phi(n, x) |d3| <= tol n^3, Phi = V(x) for even n.
n = 1: cu2qu's tangent-intersection quadratic if it passes cu2qu's test;
otherwise the s = 0 quadratic if sqrt(3)/36 |d3| <= tol and its end tangents
are within delta.

delta = 0 is "optu"; delta = infinity is "zero" (s = 0).
"""
import json, math, os
import numpy as np
from common import CORPUS, phi_max
from q1_relax import V
from schemes import optu_weights, cross
from corpus_stats import BIG, TOLS

C = np.array(json.load(open(CORPUS)), dtype=float)
N = len(C)
d3v = C[:, 3] - 3 * C[:, 2] + 3 * C[:, 1] - C[:, 0]
d3 = np.hypot(d3v[:, 0], d3v[:, 1])
ft_n = np.load(os.path.join(BIG, 'corpus_ft_n.npy'))
E0 = C[:, 1] - C[:, 0]; E1 = C[:, 2] - C[:, 3]          # handles, pointing away from the end points


def max_relax(n, tand):
    """largest 1 - x allowed by the budget at both ends (array), capped at 1."""
    out = np.ones(N)
    for E, sign in ((E0, 1.0), (E1, -1.0)):
        # start: direction E0 - (1-x) d3/(6n^2);  end: direction E1 + (1-x) d3/(6n^2)
        dv = sign * d3v
        num = 6 * n * n * tand * (E * E).sum(-1)
        den = np.abs(cross(dv, E)) + tand * (dv * E).sum(-1)
        lim = np.where(den > 1e-12, num / np.where(den > 1e-12, den, 1.0), np.inf)
        degenerate = (E * E).sum(-1) < 1e-18           # zero handle: no tangent to keep by this rule; treat as kept
        lim = np.where(degenerate, 0.0, lim)
        out = np.minimum(out, lim)
    return np.clip(out, 0.0, 1.0)


XS = np.linspace(0, 1, 21)
_odd = {}
def Phi_odd(n, x):
    """LP constant for odd n with s_0 = x (interpolated from a table over x)."""
    if n not in _odd:
        _odd[n] = np.array([phi_max(optu_weights(n, float(xx))) for xx in XS])
    return np.interp(x, XS, _odd[n])

Vv = np.vectorize(V)


def run(delta_deg, tol, j):
    n_out = np.zeros(N, dtype=int)
    left = np.ones(N, dtype=bool)
    one = ft_n[:, j] == 1
    n_out[one] = 1; left[one] = False
    if delta_deg > 0:
        tand = math.tan(math.radians(delta_deg)) if delta_deg < 89 else 1e9
        fit1 = left & (3 ** 0.5 / 36 * d3 <= tol) & (max_relax(1, tand) >= 1.0)
        n_out[fit1] = 1; left[fit1] = False
    for n in range(2, 60):
        if not left.any():
            break
        if delta_deg > 0:
            x = 1.0 - max_relax(n, tand)
        else:
            x = np.ones(N)
        if n % 2 == 0:
            Phi = Vv(x)
        else:
            Phi = Phi_odd(n, x) if n <= 15 else Vv(x) * 1.0002
        fit = left & (Phi * d3 <= tol * n ** 3)
        n_out[fit] = n; left[fit] = False
    return n_out


if __name__ == '__main__':
    print(f'{N} cubics; average quadratics per cubic, n from the equal-parameter bound; cu2qu: ' + ' / '.join(f'{ft_n[:, j].mean():.4f}' for j in range(3)))
    print('   tangent budget     tol 1.0   tol 0.5   tol 0.1     vs cu2qu')
    for delta in (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 90.0):
        means = [run(delta, tol, j).mean() for j, tol in enumerate(TOLS)]
        label = 'kept exactly' if delta == 0 else ('unbounded' if delta >= 89 else f'{delta:g} degrees')
        print(f'   {label:16s}   ' + '   '.join(f'{m:.4f}' for m in means) + '      '
              + ' '.join(f'{100 * (m / ft_n[:, j].mean() - 1):+6.1f}%' for j, m in enumerate(means)))
