#!/usr/bin/env python3
"""Look at cubics that take two quadratics by the equal-parameter test but one
by the true distance, and re-measure a few with an independent brute force."""
import numpy as np
from common import *

c = load_cubics()
g = load_gt(W + '/gt_ramp.bin', 16)
ep, hsb, hbs, fob, fos = (g[:, :, j] for j in range(5))
H = np.fmax(hsb, hbs)
tol = 1.0
n_a = first_fit(ep, tol)
n_c = first_fit(H, tol)

def brute(cub, q, M=20001):
    """two-sided Hausdorff by all-pairs on dense samples (independent of gt.c)"""
    t = np.linspace(0, 1, M)[:, None]
    B = (1-t)**3*cub[0] + 3*(1-t)**2*t*cub[1] + 3*(1-t)*t*t*cub[2] + t**3*cub[3]
    Sq = (1-t)**2*cub[0] + 2*(1-t)*t*q + t*t*cub[3]
    # chunked nearest
    def directed(A, Bc):
        best = 0
        for i in range(0, len(A), 500):
            d = np.sqrt(((A[i:i+500, None, :] - Bc[None, ::4, :])**2).sum(-1)).min(1)
            best = max(best, d.max())
        return best
    return directed(Sq, B), directed(B, Sq)

def tang_int(cub):
    a, b, cc, d = cub
    ab = b - a; cd = d - cc
    p = np.array([-ab[1], ab[0]])
    h = p.dot(a - cc) / p.dot(cd)
    return cc + cd * h

sel = np.where((n_a == 2) & (n_c == 1))[0]
print(len(sel), 'cubics with n_a = 2, n_true = 1 at tol', tol)
rng = np.random.default_rng(0)
for i in rng.choice(sel, 12, replace=False):
    q = tang_int(c[i])
    bs, bb = brute(c[i], q)
    print(i, c[i].astype(int).tolist(), f'ep1 {ep[i,0]:.3f}  gt: h_sb {hsb[i,0]:.4f} h_bs {hbs[i,0]:.4f}   brute: {bs:.4f} {bb:.4f}   fo_b {fob[i,0]:.4f} fo_s {fos[i,0]:.4f}')

# distribution of ep1 / H1 over all cubics with a tangent intersection
m = np.isfinite(ep[:, 0]) & (H[:, 0] > 1e-9)
r = ep[m, 0] / H[m, 0]
print('ep(1)/H(1) percentiles 0/5/25/50/75/95/100:', np.percentile(r, [0, 5, 25, 50, 75, 95, 100]).round(3))
print('share of cubics with H(1) <= 1:', (H[:, 0] <= 1).mean(), ' ep(1) <= 1:', (ep[:, 0] <= 1).mean())
