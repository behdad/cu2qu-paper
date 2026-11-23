#!/usr/bin/env python3
"""Q3: the tangent-intersection spline is good geometrically and bad in cu2qu's
equal-parameter measure.  On the corpus (inflection-free cubics whose pieces
all have a valid tangent intersection), at equal-parameter knots:

    equal-parameter distance / Hausdorff distance,  for T and for the cu2qu ramp.
"""
import json, os
import numpy as np
import chd, schemes
from common import CORPUS

BIG = '/home/behdad/hb/cu2qu-work/research-weights'
C = np.array(json.load(open(CORPUS)), dtype=float)
N = len(C)
convex = schemes.inflections(C) == 0
print('equal-parameter distance / Hausdorff distance (cubics with Hausdorff distance > 0.01)')
print('   scheme  n   cubics   10%     median   90%     99%')
for n in (1, 2, 3, 4):
    QT, valid = schemes.tangents(C, schemes.uniform_knots(N, n), return_valid=True)
    good = valid.all(1) & convex
    for name, Q in (('T', QT), ('ramp', schemes.ramp(C, n) if n > 1 else QT), ('optu', schemes.optu(C, n) if n > 1 else QT)):
        h = chd.hd(C, Q, m=16); e = chd.ep(C, Q, m=32)
        ok = good & (h > 0.01)
        r = e[ok] / h[ok]
        q = np.quantile(r, [0.1, 0.5, 0.9, 0.99])
        print(f'   {name:6s} {n:2d}  {ok.sum():6d}  ' + '  '.join(f'{v:7.2f}' for v in q))
