#!/usr/bin/env python3
"""Where does the closed form give away quadratics?  For cubics where it picks a
larger n than the truth: is the bound invalid at the true n, or just loose?"""
import sys
import numpy as np
from common import *
from closed import *

c = load_cubics()
NMAX = 16
g = load_gt(W + '/gt_ramp.bin', NMAX)
H = np.fmax(g[:, :, 1], g[:, :, 2])
ep = g[:, :, 0]
N = len(c)
B = np.full((N, NMAX), np.inf); Braw = np.full((N, NMAX), np.inf); V = np.zeros((N, NMAX), bool)
b1, v1, det = n1_closed(c, detail=True)
Braw[:, 0] = b1; V[:, 0] = v1
B[:, 0] = np.where(v1, np.minimum(b1, ep[:, 0]), ep[:, 0])
d3 = np.hypot(*d3_of(c).T)
parts = {}
for n in range(2, NMAX + 1):
    bn, vn, pp = ramp_closed(c, n, detail=True)
    parts[n] = pp
    epn = d3 * Phi_ramp(n) / n ** 3
    Braw[:, n - 1] = bn; V[:, n - 1] = vn
    B[:, n - 1] = np.where(vn, np.minimum(bn, epn), epn)

for tol in (1.0, 0.5, 0.1):
    na = first_fit(ep, tol); nb = first_fit(B, tol); nc = first_fit(H, tol * (1 + 1e-9))
    lost = nb > nc
    print(f'tol {tol}: cubics where closed form takes more than the truth: {lost.sum()} ({100*lost.mean():.2f}%), extra quadratics {int((nb-nc)[lost].sum())}')
    for n in range(1, 7):
        sel = lost & (nc == n)
        if not sel.any():
            continue
        inv = sel & ~V[:, n - 1]
        loose = sel & V[:, n - 1]
        r = Braw[loose, n - 1] / H[loose, n - 1]
        txt = f'   true n = {n}: {sel.sum():5d} cubics;  invalid at n: {inv.sum():5d};  valid but loose: {loose.sum():5d}'
        if loose.any():
            txt += f'  (bound/H: median {np.median(r):.3f}, 90% {np.percentile(r, 90):.3f}; bound/tol median {np.median(Braw[loose, n-1])/tol:.3f})'
        print(txt)
        if n == 1 and inv.any():
            l1, l2 = det['l1'], det['l2']
            c_l = inv & ~((l1 > 0) & (l1 < 2) & (l2 > 0) & (l2 < 2))
            c_ab = inv & ~((det['alp'] > 0) & (det['bep'] > 0))
            c_x = inv & ~(det['xmax'] <= 0.75)
            c_arc = inv & ~det['arc']
            print(f'        n=1 invalid because: handles outside (0,2): {c_l.sum()}, minorant <= 0: {c_ab.sum()}, discriminant: {c_x.sum()}, arc: {c_arc.sum()}')
        if n >= 2 and loose.any():
            pp = parts[n][loose]
            k = np.argmax(pp.sum(-1), axis=1)
            T = pp[np.arange(len(k)), k]
            print(f'        dominating piece: first/last {np.mean((k == 0) | (k == n - 1)):.2f};  T2 share of bound: median {np.median(T[:, 1] / T.sum(1)):.3f}, 90% {np.percentile(T[:, 1] / T.sum(1), 90):.3f};  T3 share median {np.median(T[:, 2] / T.sum(1)):.4f}')
