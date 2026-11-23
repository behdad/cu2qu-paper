#!/usr/bin/env python3
"""Is the yardstick optimiser (opt.py) finding the optimum?  An independent
check where the search space has two dimensions and can be scanned:

  n = 2, tangents kept:  Q_0 = c0 + a e0, Q_1 = c3 - b e1   (a, b)
  n = 1, tangents free:  Q anywhere                          (x, y)

A coarse global grid, then four rounds of local grid refinement, minimising the
two-sided Hausdorff distance directly (chd.hd1); compared with the optimiser's
value stored by yardstick_run.py.
"""
import json, os
import numpy as np
import chd
from opt import end_tangents
from common import CORPUS
from canon import SHAPES

BIG = '/home/behdad/hb/cu2qu-work/research-weights'


def refine(f, x0, y0, sx, sy, rounds=5, k=10):
    best = (f(x0, y0), x0, y0)
    for r in range(rounds):
        for i in range(-k, k + 1):
            for j in range(-k, k + 1):
                x = best[1] + i * sx / k if r else x0 + i * sx / k
                y = best[2] + j * sy / k if r else y0 + j * sy / k
                v = f(x, y)
                if v < best[0]:
                    best = (v, x, y)
        sx /= 4; sy /= 4
        x0, y0 = best[1], best[2]
    return best


def tan2(c):
    e0, e1 = end_tangents(c)
    L = np.hypot(*(c[3] - c[0]))
    def f(a, b):
        if a <= 0 or b <= 0:
            return 1e30
        return chd.hd1(c, np.array([c[0] + a * e0, c[3] - b * e1]), m=24)[0]
    # global coarse scan
    best = (1e30, 0, 0)
    for a in np.linspace(0.02, 1.5, 40) * L:
        for b in np.linspace(0.02, 1.5, 40) * L:
            v = f(a, b)
            if v < best[0]:
                best = (v, a, b)
    return refine(f, best[1], best[2], 0.04 * L, 0.04 * L)[0]


def free1(c):
    lo = c.min(0); hi = c.max(0); span = (hi - lo).max() + 1e-9
    def f(x, y):
        return chd.hd1(c, np.array([[x, y]]), m=32)[0]
    best = (1e30, 0, 0)
    for x in np.linspace(lo[0] - 0.3 * span, hi[0] + 0.3 * span, 60):
        for y in np.linspace(lo[1] - 0.3 * span, hi[1] + 0.3 * span, 60):
            v = f(x, y)
            if v < best[0]:
                best = (v, x, y)
    return refine(f, best[1], best[2], 0.03 * span, 0.03 * span)[0]


if __name__ == '__main__':
    C = np.array(json.load(open(CORPUS)), dtype=float)
    T = np.load(os.path.join(BIG, 'corpus_yard_tan.npy'))
    Fr = np.load(os.path.join(BIG, 'corpus_yard_free.npy'))
    rng = np.random.default_rng(8)
    ok = np.where(np.isfinite(T[:, 1]) & (T[:, 1] > 0.05))[0]
    pick = rng.choice(ok, 40, replace=False)
    r2 = []; r1 = []
    for i in pick:
        g2 = tan2(C[i]); r2.append(T[i, 1] / g2)
    ok1 = np.where(np.isfinite(Fr[:, 0]) & (Fr[:, 0] > 0.05))[0]
    pick1 = rng.choice(ok1, 40, replace=False)
    for i in pick1:
        g1 = free1(C[i]); r1.append(Fr[i, 0] / g1)
    r2 = np.array(r2); r1 = np.array(r1)
    print('optimiser value / grid-search value (1 = same; above 1 = the optimiser missed something)')
    print(f'  n = 2, tangents kept, 40 corpus cubics: min {r2.min():.4f}  median {np.median(r2):.4f}  max {r2.max():.4f}')
    print(f'  n = 1, tangents free, 40 corpus cubics: min {r1.min():.4f}  median {np.median(r1):.4f}  max {r1.max():.4f}')
    from opt import best_spline
    import schemes
    print('  canonical shapes (grid, optimiser):')
    for name, pts in SHAPES.items():
        c = np.asarray(pts, dtype=float); Cc = c[None]
        g2 = tan2(c)
        st = [schemes.ramp(Cc, 2)[0], schemes.tangents(Cc, schemes.uniform_knots(1, 2))[0], schemes.tangents(Cc, schemes.affine_knots(Cc, 2))[0]]
        o2 = best_spline(c, 2, 'tan', st)[0]
        g1 = free1(c)
        o1 = best_spline(c, 1, 'free', [schemes.zero(Cc, 1)[0], schemes.tangents(Cc, schemes.uniform_knots(1, 1))[0]])[0]
        print(f'    {name:28s} n=2 kept: {g2:9.4f} {o2:9.4f}   n=1 free: {g1:9.4f} {o1:9.4f}')
