#!/usr/bin/env python3
"""Q4: how far the closed-form geometric constructions are from the numerical
optimum, at equal n, on the whole corpus (inflection-free cubics whose pieces
all have a valid tangent intersection; yardstick error above 0.02 units).

  T        tangent-intersection spline (tangents kept)            vs TAN*
  BT-kept  balanced T, tangents kept                               vs TAN*
  BT-free  balanced T, tangents free                               vs FREE*
  single   one quadratic, Q_T - 2 (sqrt 2 - 1) P_N (Q_T - q_M)     vs FREE*(1)
  zero     d3 family, s = 0                                        vs FREE*
at gamma knots (closed form) and at equal affine-arc knots (numerical).
"""
import json, os
import numpy as np
import chd, schemes, balanced
from common import CORPUS

BIG = '/home/behdad/hb/cu2qu-work/research-weights'
C = np.array(json.load(open(CORPUS)), dtype=float)
N = len(C)
T = np.load(os.path.join(BIG, 'corpus_yard_tan.npy'))
Fr = np.load(os.path.join(BIG, 'corpus_yard_free.npy'))
convex = schemes.inflections(C) == 0


def line(label, H, yard, n, mask):
    ok = mask & np.isfinite(H) & np.isfinite(yard[:, n - 1]) & (yard[:, n - 1] > 0.02)
    if ok.sum() < 50:
        return f'{label:28s} n={n}:      -'
    r = H[ok] / yard[ok, n - 1]
    q = np.quantile(r, [0.1, 0.5, 0.9, 0.99])
    return f'{label:28s} n={n}: cubics {ok.sum():6d}   ratio to the optimum: 10% {q[0]:5.2f}  median {q[1]:5.2f}  90% {q[2]:5.2f}  99% {q[3]:6.2f}'


if __name__ == '__main__':
    for kn_label, knots_of in (('gamma knots', lambda n: schemes.gamma_knots(C, n)), ('affine-arc knots', lambda n: schemes.affine_knots(C, n))):
        print('==', kn_label)
        for n in (1, 2, 3, 4):
            kn = knots_of(n)
            QT, valid = schemes.tangents(C, kn, return_valid=True)
            good = valid.all(1) & convex
            print(line('T', chd.hd(C, QT, m=16), T, n, good))
            Qk, vk, _ = balanced.bt(C, kn, True)
            print(line('BT-kept', chd.hd(C, Qk, m=16), T, n, good & vk))
            Qf, vf, _ = balanced.bt(C, kn, False)
            print(line('BT-free', chd.hd(C, Qf, m=16), Fr, n, good & vf))
            if n == 1:
                Qs, vs, _ = balanced.single(C, 'N')
                print(line('single (normal rule)', chd.hd(C, Qs, m=16), Fr, n, good & vs))
            print(line('zero (equal-parameter)', chd.hd(C, schemes.zero(C, n), m=16), Fr, n, good))
        print()
