#!/usr/bin/env python3
"""Q4 on the corpus: for every cubic, the numerically best implied-on-curve
spline for n = 1, 2, ... (stopping once below STOP), with both end tangents
kept (TAN*) and with tangents free (FREE*).

Writes BIG/corpus_yard_tan.npy and BIG/corpus_yard_free.npy, (N, NMAX)
Hausdorff distances, NaN where not computed.

The values are upper bounds on the true optimum (a local optimiser with
several starts); every closed-form scheme is among the starts, so the
yardstick is never worse than any of them.

Usage: yardstick_run.py [count [seed]]     (default: the whole corpus)
"""
import json, os, sys, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from multiprocessing import Pool
import chd, schemes
from opt import best_spline
from common import CORPUS

BIG = '/home/behdad/hb/cu2qu-work/research-weights'
NMAX = 24
NY = 10               # largest n optimised
STOP = 0.1


def one(c):
    C = c[None]
    tan = np.full(NMAX, np.nan); free = np.full(NMAX, np.nan)
    need_t = True; need_f = True
    for n in range(1, NY + 1):
        if not (need_t or need_f):
            break
        Tu, vu = schemes.tangents(C, schemes.uniform_knots(1, n), return_valid=True)
        Ta, va = schemes.tangents(C, schemes.affine_knots(C, n), return_valid=True)
        Tg, vg = schemes.tangents(C, schemes.gamma_knots(C, n), return_valid=True)
        Qr = schemes.ramp(C, n)[0] if n > 1 else None
        Qo = schemes.optu(C, n)[0] if n > 1 else None
        Qz = schemes.zero(C, n)[0]
        qt = None
        if need_t:
            if n == 1:
                if vu.all():
                    tan[0] = chd.hd1(c, Tu[0], m=48)[0]; qt = Tu[0]
                else:
                    tan[0] = np.inf
            else:
                starts = [Qo, Qr, Tu[0], Ta[0], Tg[0]]      # T with fallbacks still keeps the end tangents
                r = best_spline(c, n, 'tan', starts)
                tan[n - 1] = r[0]; qt = r[1]
            if tan[n - 1] <= STOP:
                need_t = False
        if need_f:
            starts = [Qz] + ([qt] if qt is not None else []) + ([Ta[0]] if va.all() else []) + ([Qo] if Qo is not None else [])
            r = best_spline(c, n, 'free', starts)
            free[n - 1] = r[0]
            if need_t is False and qt is not None:
                free[n - 1] = min(free[n - 1], tan[n - 1])
            if free[n - 1] <= STOP:
                need_f = False
    return tan, free


def work(args):
    i0, Cs = args
    out_t = np.empty((len(Cs), NMAX)); out_f = np.empty((len(Cs), NMAX))
    for i, c in enumerate(Cs):
        try:
            out_t[i], out_f[i] = one(c)
        except Exception as ex:           # keep going; report
            out_t[i] = np.nan; out_f[i] = np.nan
            print('failed', i0 + i, c.tolist(), repr(ex), flush=True)
    return i0, out_t, out_f


if __name__ == '__main__':
    C = np.array(json.load(open(CORPUS)), dtype=float)
    idx = np.arange(len(C))
    tag = ''
    if len(sys.argv) > 1:
        count = int(sys.argv[1]); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
        idx = np.sort(np.random.default_rng(seed).choice(len(C), count, replace=False))
        tag = f'_sample{count}'
    Cs = C[idx]
    jobs = [(i, Cs[i:i + 200]) for i in range(0, len(Cs), 200)]
    T = np.full((len(Cs), NMAX), np.nan); Fr = np.full((len(Cs), NMAX), np.nan)
    t0 = time.time(); done = 0
    with Pool(int(os.environ.get("YARD_PROCS", "22"))) as pool:
        for i0, a, b in pool.imap_unordered(work, jobs):
            T[i0:i0 + len(a)] = a; Fr[i0:i0 + len(b)] = b
            done += len(a)
            if done % 4000 < 200:
                print(f'{done} / {len(Cs)}  ({time.time() - t0:.0f} s)', flush=True)
    np.save(os.path.join(BIG, f'corpus_yard_tan{tag}.npy'), T)
    np.save(os.path.join(BIG, f'corpus_yard_free{tag}.npy'), Fr)
    if tag:
        np.save(os.path.join(BIG, f'corpus_yard_idx{tag}.npy'), idx)
    print('done', time.time() - t0)
