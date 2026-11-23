#!/usr/bin/env python3
"""Geometric (Hausdorff) error of each closed-form scheme on every corpus cubic,
for n = 1, 2, ... until the error is below the smallest tolerance of interest.

Writes  BIG/corpus_<scheme>.npy : (N, NMAX) Hausdorff distances (NaN = not
computed because a smaller n already got below STOP), and for the d3-family
schemes also the equal-parameter error.

Usage: corpus_run.py [scheme ...]
"""
import json, os, sys, time
import numpy as np
import chd, schemes
from common import CORPUS

BIG = '/home/behdad/hb/cu2qu-work/research-weights'
NMAX = 24
STOP = 0.05          # stop refining a cubic once below this (tolerances studied: 1, 0.5, 0.1)
M = 16               # samples per quadratic for the instrument


def scheme_points(name, C, n):
    N = len(C)
    if name == 'ramp':      # cu2qu; n = 1 is the tangent intersection
        return schemes.ramp(C, n) if n > 1 else schemes.tangents(C, schemes.uniform_knots(N, 1))
    if name == 'optu':
        return schemes.optu(C, n) if n > 1 else schemes.tangents(C, schemes.uniform_knots(N, 1))
    if name == 'zero':
        return schemes.zero(C, n)
    if name.startswith('relax'):
        x = int(name[5:]) / 100.0
        return schemes.relaxed(C, n, x) if n > 1 else schemes.zero(C, 1)
    if name == 'Tuni':
        return schemes.tangents(C, schemes.uniform_knots(N, n))
    if name == 'Taff':
        return schemes.tangents(C, schemes.affine_knots(C, n))
    if name == 'Tgam':
        return schemes.tangents(C, schemes.gamma_knots(C, n))
    raise ValueError(name)


def run(name, C):
    N = len(C)
    H = np.full((N, NMAX), np.nan)
    active = np.arange(N)
    t0 = time.time()
    for n in range(1, NMAX + 1):
        if len(active) == 0:
            break
        Q = scheme_points(name, C[active], n)
        h = chd.hd(C[active], Q, m=M)
        H[active, n - 1] = h
        active = active[~(h <= STOP)]
        print(f'  {name} n = {n:2d}: {N - len(active):6d} done, {len(active):6d} left  ({time.time() - t0:.0f} s)', flush=True)
    return H


if __name__ == '__main__':
    C = np.array(json.load(open(CORPUS)), dtype=float)
    names = sys.argv[1:] or ['ramp', 'optu', 'zero', 'relax50', 'Tuni', 'Taff', 'Tgam']
    for name in names:
        H = run(name, C)
        np.save(os.path.join(BIG, f'corpus_{name}.npy'), H)
