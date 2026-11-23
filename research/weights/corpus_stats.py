#!/usr/bin/env python3
"""Average quadratics per cubic on the corpus, for each scheme, at tolerance
1.0, 0.5 and 0.1 font units.

Two ways of choosing n:
  formula  the least n allowed by the scheme's closed-form equal-parameter bound
           Phi(n) |d3| <= tol n^3 (d3-family schemes only); n = 1 is tried first
           the way cu2qu does (tangent intersection, equal-parameter test)
  oracle   the least n whose spline is within tol of the cubic in Hausdorff
           distance, measured (corpus_run.py)

Needs corpus_run.py to have been run.  Yardstick rows (TAN*, FREE*) appear if
yardstick_run.py has been run.
"""
import json, os, sys
import numpy as np
from common import CORPUS, FT_LIB, phi_max, ramp as ramp_w
from schemes import optu_weights

BIG = '/home/behdad/hb/cu2qu-work/research-weights'
TOLS = (1.0, 0.5, 0.1)


def least_n(H, tol):
    """H (N, NMAX) -> least n (1-based) with H <= tol; NMAX + 1 where none."""
    ok = H <= tol
    first = ok.argmax(1) + 1
    first[~ok.any(1)] = H.shape[1] + 1
    return first


def formula_n(d3, tol, Phi, nmax=200, start=1):
    """least n >= start with Phi(n) d3 <= tol n^3."""
    out = np.zeros(len(d3), dtype=int)
    left = np.ones(len(d3), dtype=bool)
    for n in range(start, nmax + 1):
        fit = left & (Phi(n) * d3 <= tol * n ** 3)
        out[fit] = n
        left &= ~fit
        if not left.any():
            break
    return out


if __name__ == '__main__':
    C = np.array(json.load(open(CORPUS)), dtype=float)
    N = len(C)
    d3v = C[:, 3] - 3 * C[:, 2] + 3 * C[:, 1] - C[:, 0]
    d3 = np.hypot(d3v[:, 0], d3v[:, 1])

    cache = os.path.join(BIG, 'corpus_ft_n.npy')
    if os.path.exists(cache):
        ft_n = np.load(cache)
    else:
        sys.path.insert(0, FT_LIB)
        from fontTools.cu2qu.cu2qu import curve_to_quadratic
        ft_n = np.zeros((N, len(TOLS)), dtype=int)
        for i in range(N):
            cc = [tuple(p) for p in C[i]]
            for j, tol in enumerate(TOLS):
                ft_n[i, j] = len(curve_to_quadratic(cc, tol)) - 2
        np.save(cache, ft_n)

    Phi_ramp = lambda n: phi_max(ramp_w(n))
    _opt = {}
    def Phi_opt(n):
        if n not in _opt:
            _opt[n] = 4 / 27 if n % 2 == 0 else (phi_max(optu_weights(n)) if n <= 31 else 4 / 27 * (1 + 1e-9))
        return _opt[n]
    def Phi_even(n):          # alternating on even n, (1, 0, -1) for n = 3, no other odd n
        return 4 / 27 if n % 2 == 0 else (27 / 128 if n == 3 else float('inf'))
    Phi_zero = lambda n: 3 ** 0.5 / 36
    _r50 = {}
    def Phi_r50(n):
        if n not in _r50:
            _r50[n] = phi_max(optu_weights(n, 0.5)) if n > 1 else 3 ** 0.5 / 36
        return _r50[n]

    rows = []
    def add(label, how, ns):
        rows.append((label, how, [float(np.mean(ns[j])) for j in range(len(TOLS))], ns))

    add('cu2qu (fontTools itself)', 'as implemented', [ft_n[:, j] for j in range(len(TOLS))])
    one = [ft_n[:, j] == 1 for j in range(len(TOLS))]
    def with_single(ns, j):
        return np.where(one[j], 1, ns)
    f_ramp = [with_single(formula_n(d3, tol, Phi_ramp, start=2), j) for j, tol in enumerate(TOLS)]
    # the formula and fontTools' own trial agree except on exact ties (|d3| Phi(n) = tol n^3
    # in exact arithmetic, e.g. |d3| = 128 at tolerance 1), where fontTools' recursive test rounds the other way
    for j, tol in enumerate(TOLS):
        bad = np.where(f_ramp[j] != ft_n[:, j])[0]
        for b in bad:
            n = f_ramp[j][b]
            assert abs(Phi_ramp(n) * d3[b] - tol * n ** 3) < 1e-9 * tol * n ** 3, (C[b], n, ft_n[b, j])
        print(f'tolerance {tol}: closed-form n for the ramp differs from fontTools on {len(bad)} cubics, all exact ties')
    add('ramp', 'formula', f_ramp)
    add('optu (alternating / best odd)', 'formula', [with_single(formula_n(d3, tol, Phi_opt, start=2), j) for j, tol in enumerate(TOLS)])
    def Phi_block(n):         # pairs (+1, -1), plus one block of 3 (n = 3) or of 5 (odd n >= 5)
        return 4 / 27 if n % 2 == 0 else (27 / 128 if n == 3 else 0.151879364835661353)
    add('pairs + one block of 5 for odd n', 'formula', [with_single(formula_n(d3, tol, Phi_block, start=2), j) for j, tol in enumerate(TOLS)])
    add('alternating, even n and 3 only', 'formula', [with_single(formula_n(d3, tol, Phi_even, start=2), j) for j, tol in enumerate(TOLS)])
    add('relax50 (s = +-0.5 at the ends)', 'formula', [formula_n(d3, tol, Phi_r50) for tol in TOLS])
    add('zero (s = 0)', 'formula', [formula_n(d3, tol, Phi_zero) for tol in TOLS])

    names = [('ramp', 'ramp'), ('optu', 'optu (alternating / best odd)'), ('Tuni', 'T at equal-parameter knots'),
             ('Tgam', 'T at closed-form (gamma) knots'), ('Taff', 'T at equal affine-arc knots'),
             ('relax50', 'relax50'), ('zero', 'zero (s = 0)'),
             ('yard_tan', 'TAN* (numerical optimum, tangents kept)'), ('yard_free', 'FREE* (numerical optimum, tangents free)')]
    Hs = {}
    for key, label in names:
        path = os.path.join(BIG, f'corpus_{key}.npy')
        if not os.path.exists(path):
            continue
        H = np.load(path)
        Hs[key] = H
        add(label, 'oracle', [least_n(H, tol) for tol in TOLS])

    print(f'{N} cubics.  Average quadratics per cubic:')
    print(f'  {"scheme":44s} {"n chosen by":15s} ' + ''.join(f'tol {t:<8}' for t in TOLS) + '   vs cu2qu')
    base = rows[0][2]
    for label, how, means, ns in rows:
        print(f'  {label:44s} {how:15s} ' + ''.join(f'{m:<12.4f}' for m in means)
              + '   ' + ' '.join(f'{100 * (m / b - 1):+6.1f}%' for m, b in zip(means, base)))
    print()
    print('Distribution of n (percent of cubics), tolerance 1.0 / 0.5 / 0.1:')
    for label, how, means, ns in rows:
        txt = []
        for j in range(len(TOLS)):
            cnt = np.bincount(np.minimum(ns[j], 9), minlength=10)
            txt.append(' '.join(f'{100 * cnt[k] / N:5.1f}' for k in range(1, 10)))
        print(f'  {label:44s} {how:15s}')
        for j, tol in enumerate(TOLS):
            print(f'      tol {tol:<4}:  n=1..8, 9+:  {txt[j]}')
