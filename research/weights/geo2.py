#!/usr/bin/env python3
"""Q5: every scheme with n chosen by the same closed-form geometric ESTIMATE
(gest.g_estimate, 17 nodes per piece), then checked against the measured
Hausdorff distance.

    n = least n <= NTRY with SAFETY * estimate(n) <= tol;
    otherwise the scheme's d3-family fallback with its equal-parameter formula.

Schemes
  kept:  ramp, optu, T (tangent intersections at gamma knots; d3-family point in
         a piece without a valid intersection), BT-kept (balanced T; plain T
         where a piece is not valid)
  free:  zero, BT-free (balanced T; zero where not valid)

Output: average quadratics per cubic; share of cubics whose chosen spline is in
fact farther than tol; median and worst excess.
"""
import json, os, sys
import numpy as np
import chd, schemes, balanced
from gest import g_estimate
from common import CORPUS, phi_max, ramp as ramp_w
from schemes import optu_weights
from corpus_stats import formula_n, BIG, TOLS

SAFETY = float(os.environ.get('SAFETY', '1.1'))
NTRY = 10
NODES = 17

C = np.array(json.load(open(CORPUS)), dtype=float)
N = len(C)
d3v = C[:, 3] - 3 * C[:, 2] + 3 * C[:, 1] - C[:, 0]
d3 = np.hypot(d3v[:, 0], d3v[:, 1])
ft_n = np.load(os.path.join(BIG, 'corpus_ft_n.npy'))
Phi_ramp = lambda n: phi_max(ramp_w(n))
Phi_opt = lambda n: 4 / 27 if n % 2 == 0 else phi_max(optu_weights(n))
Phi_zero = lambda n: 3 ** 0.5 / 36


def build(kind, Cs, n):
    """(off-curve points, knots used, usable).  Tangent intersections are used only for
    cubics without an inflection whose pieces all have a valid intersection; the
    rest get the d3-family points (optu if tangents are kept, s = 0 if not).  With
    tangents kept there is no single quadratic for such a cubic (usable = False)."""
    M = len(Cs)
    uk = schemes.uniform_knots(M, n)
    usable = np.ones(M, dtype=bool)
    if kind == 'zero':
        return schemes.zero(Cs, n), uk, usable
    convex = schemes.inflections(Cs) == 0
    if kind in ('ramp', 'optu'):
        if n == 1:
            Q, valid = schemes.tangents(Cs, uk, return_valid=True)
            return Q, uk, valid.all(1) & convex
        return (schemes.ramp(Cs, n) if kind == 'ramp' else schemes.optu(Cs, n)), uk, usable
    kn = schemes.gamma_knots(Cs, n)
    QT, valid = schemes.tangents(Cs, kn, return_valid=True)
    good = valid.all(1) & convex
    if kind in ('T', 'BT-kept'):
        Q = QT
        if kind == 'BT-kept':
            Qb, vb, est = balanced.bt(Cs, kn, True)
            Q = np.where((good & vb)[:, None, None], Qb, QT)
        if n == 1:
            return Q, kn, good
        fb = schemes.optu(Cs, n)
        return np.where(good[:, None, None], Q, fb), np.where(good[:, None], kn, uk), usable
    if kind == 'BT-free':
        Qb, vb, est = balanced.bt(Cs, kn, False)
        good = good & vb
        fb = schemes.zero(Cs, n)
        return np.where(good[:, None, None], Qb, fb), np.where(good[:, None], kn, uk), usable
    raise ValueError(kind)


def run(kind, tol, jtol, want_groups=False):
    kept = kind in ('ramp', 'optu', 'T', 'BT-kept')
    if kept:
        Phi = Phi_ramp if kind == 'ramp' else Phi_opt
        n_fb = np.where(ft_n[:, jtol] == 1, 1, formula_n(d3, tol, Phi, start=2))
    else:
        n_fb = formula_n(d3, tol, Phi_zero)
    n_out = n_fb.copy()
    H = np.full(N, np.nan)
    left = np.ones(N, dtype=bool)
    groups = []
    for n in range(1, NTRY + 1):
        cand = np.where(left & (n_fb > n))[0]          # the estimate can only save quadratics over the formula
        if len(cand) == 0:
            break
        Q, kn, usable = build(kind, C[cand], n)
        est = g_estimate(C[cand], Q, kn, NODES)
        okk = usable & np.isfinite(est) & (SAFETY * est <= tol)
        take = cand[okk]
        n_out[take] = n; left[take] = False
        if len(take):
            H[take] = chd.hd(C[take], Q[okk], m=16)
            groups.append((take, Q[okk]))
    by_est = ~left
    rest = np.where(left)[0]
    fbk = {'ramp': 'ramp', 'optu': 'optu', 'T': 'optu', 'BT-kept': 'optu', 'zero': 'zero', 'BT-free': 'zero'}[kind]
    for n in np.unique(n_out[rest]):
        sel = rest[n_out[rest] == n]
        if n == 1 and fbk != 'zero':          # cu2qu's own single quadratic passed its test
            Q = schemes.tangents(C[sel], schemes.uniform_knots(len(sel), 1))
        else:
            Q = build(fbk, C[sel], int(n))[0]
        H[sel] = chd.hd(C[sel], Q, m=16)
        groups.append((sel, Q))
    if want_groups:
        return n_out, by_est, H, groups
    return n_out, by_est, H


if __name__ == '__main__':
    print(f'{N} cubics; safety factor on the estimate: {SAFETY}; cu2qu average n: ' + ' / '.join(f'{ft_n[:, j].mean():.4f}' for j in range(3)))
    print(f'  {"scheme":9s} {"tol":>5s}  avg n   vs cu2qu  chosen by estimate | over tolerance: cubics  share     median H/tol  worst H/tol | percent with n = 1, 2, 3, 4, 5+')
    for kind in ('ramp', 'optu', 'T', 'BT-kept', 'zero', 'BT-free'):
        for j, tol in enumerate(TOLS):
            n_out, by_est, H = run(kind, tol, j)
            over = H > tol * (1 + 1e-9)
            ex = H[over] / tol
            cnt = np.bincount(np.minimum(n_out, 5), minlength=6)[1:] * 100.0 / N
            print(f'  {kind:9s} {tol:5.1f}  {n_out.mean():.4f}  {100 * (n_out.mean() / ft_n[:, j].mean() - 1):+6.1f}%   {by_est.mean():.3f}             |'
                  f'   {over.sum():6d}  {over.mean():.5f}   {np.median(ex) if len(ex) else 1.0:8.3f}     {ex.max() if len(ex) else 1.0:8.3f}  | '
                  + ' '.join(f'{v:5.1f}' for v in cnt), flush=True)
        print()
