#!/usr/bin/env python3
"""Q3 on the corpus, in cu2qu's own (equal-parameter) measure: what the best
weights TOGETHER with sliding the end points along the end tangents would buy.

The constant Phi_slide(n; psi_0, psi_1) of q3_slide_ep.py depends on the cubic
only through the two angles from d3 to the end tangents.  It is tabulated on a
10-degree grid and interpolated (so the n below are approximate, not exact),
then  n = least n >= 2 with Phi_slide |d3| <= tol n^3;  n = 1 as cu2qu does.
"""
import json, math, os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from multiprocessing import Pool
from common import CORPUS
from schemes import unit_end_tangents, cross
from corpus_stats import BIG, TOLS

P0 = np.arange(0, 91, 10)                 # psi_0 >= 0 (mirror symmetry)
P1 = np.arange(-90, 91, 10)               # psi_1; -90 and 90 are the same line
NS = list(range(2, 11))


def job(args):
    from q3_slide_ep import solve
    n, a0, a1 = args
    t0 = math.inf if abs(a0) == 90 else math.tan(math.radians(a0))
    t1 = math.inf if abs(a1) == 90 else math.tan(math.radians(a1))
    return n, a0, a1, solve(n, t0, t1, grid=17, ndir=96)[0]


if __name__ == '__main__':
    cache = os.path.join(BIG, 'slide_table.npy')
    if os.path.exists(cache):
        tab = np.load(cache)
    else:
        jobs = [(n, int(a0), int(a1)) for n in NS for a0 in P0 for a1 in P1]
        tab = np.zeros((len(NS), len(P0), len(P1)))
        with Pool(12) as pool:
            for n, a0, a1, v in pool.imap_unordered(job, jobs, chunksize=4):
                tab[NS.index(n), list(P0).index(a0), list(P1).index(a1)] = v
        np.save(cache, tab)
    C = np.array(json.load(open(CORPUS)), dtype=float)
    N = len(C)
    d3v = C[:, 3] - 3 * C[:, 2] + 3 * C[:, 1] - C[:, 0]
    d3 = np.hypot(d3v[:, 0], d3v[:, 1])
    ft_n = np.load(os.path.join(BIG, 'corpus_ft_n.npy'))
    e0, e1 = unit_end_tangents(C)
    dh = d3v / np.where(d3 > 1e-12, d3, 1.0)[:, None]
    def angle(e):
        a = np.degrees(np.arctan2(cross(dh, e), (dh * e).sum(-1)))        # (-180, 180]
        a = np.where(a > 90, a - 180, np.where(a <= -90, a + 180, a))     # a line: (-90, 90]
        return a
    a0 = angle(e0); a1 = angle(e1)
    flip = a0 < 0
    a0 = np.where(flip, -a0, a0); a1 = np.where(flip, -a1, a1)
    a1 = np.where(a1 <= -90, a1 + 180, a1)

    def lookup(ni):
        i = np.clip(a0 / 10, 0, len(P0) - 1 - 1e-9); j = np.clip((a1 + 90) / 10, 0, len(P1) - 1 - 1e-9)
        i0 = np.floor(i).astype(int); j0 = np.floor(j).astype(int)
        fi = i - i0; fj = j - j0
        T = tab[ni]
        return (T[i0, j0] * (1 - fi) * (1 - fj) + T[i0 + 1, j0] * fi * (1 - fj)
                + T[i0, j0 + 1] * (1 - fi) * fj + T[i0 + 1, j0 + 1] * fi * fj)

    print('Phi_slide at a few angles (psi_0, psi_1):')
    for (x, y) in ((0, 0), (30, -30), (50, -50), (50, 50), (70, -20), (90, 90)):
        print(f'   ({x:3d}, {y:4d}):  ' + '  '.join(f'n={n}: {tab[NS.index(n), list(P0).index(x), list(P1).index(y)]:.4f}' for n in (2, 3, 4, 5, 6)))
    print()
    print(f'{N} cubics; average quadratics per cubic with the equal-parameter measure, tangents kept:')
    print('   cu2qu                              ' + '   '.join(f'{ft_n[:, j].mean():.4f}' for j in range(3)))
    for j, tol in enumerate(TOLS):
        pass
    means = []
    for j, tol in enumerate(TOLS):
        n_out = np.where(ft_n[:, j] == 1, 1, 0)
        left = n_out == 0
        for ni, n in enumerate(NS):
            Phi = lookup(ni)
            fit = left & (Phi * d3 <= tol * n ** 3)
            n_out[fit] = n; left &= ~fit
        # beyond the table: the d3-family optimum (sliding helps less and less? no: use 4/27 as a safe value)
        for n in range(NS[-1] + 1, 80):
            fit = left & (4 / 27 * 1.03 * d3 <= tol * n ** 3)
            n_out[fit] = n; left &= ~fit
        means.append(n_out.mean())
    print('   best weights + sliding (approx.)   ' + '   '.join(f'{m:.4f}' for m in means)
          + '     ' + ' '.join(f'{100 * (m / ft_n[:, j].mean() - 1):+6.1f}%' for j, m in enumerate(means)))
