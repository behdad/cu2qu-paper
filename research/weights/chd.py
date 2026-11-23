#!/usr/bin/env python3
"""ctypes wrapper for chd.c (Hausdorff distance, cubic against quadratic spline).

    hd(C, Q, m=24, both=True)     C (N, 4, 2), Q (N, n, 2) -> (N,)
    hd1(c, Q, m=48)               one cubic -> (hausdorff, spline->cubic, cubic->spline, where)
    ep(C, Q, m=32)                equal-parameter distance, (N,)
    resid(c, Q, m)                signed distances + gradients at xi = i/m

Running this file checks the C code against the numpy reference in geom.py and
against a dense polyline brute force.
"""
import ctypes, os, subprocess
import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_so = os.path.join(_here, 'libchd.so')
if not os.path.exists(_so) or os.path.getmtime(_so) < os.path.getmtime(os.path.join(_here, 'chd.c')):
    subprocess.check_call(['gcc', '-O2', '-fopenmp', '-shared', '-fPIC', '-o', _so,
                           os.path.join(_here, 'chd.c'), '-lm'])
_lib = ctypes.CDLL(_so)
_dp = ctypes.POINTER(ctypes.c_double)
_lib.hd_one.restype = ctypes.c_double
_lib.hd_one.argtypes = [_dp, _dp, ctypes.c_int, ctypes.c_int, ctypes.c_int, _dp, _dp, _dp]
_lib.hd_batch.restype = None
_lib.hd_batch.argtypes = [_dp, _dp, ctypes.c_long, ctypes.c_int, ctypes.c_int, ctypes.c_int, _dp]
_lib.ep_batch.restype = None
_lib.ep_batch.argtypes = [_dp, _dp, ctypes.c_long, ctypes.c_int, ctypes.c_int, _dp]
_lib.resid.restype = None
_lib.resid.argtypes = [_dp, _dp, ctypes.c_int, ctypes.c_int, _dp, _dp, _dp, _dp]


def _p(a):
    return a.ctypes.data_as(_dp)


def hd(C, Q, m=24, both=True):
    C = np.ascontiguousarray(C, dtype=np.float64).reshape(-1, 8)
    N = C.shape[0]
    Q = np.ascontiguousarray(Q, dtype=np.float64).reshape(N, -1)
    n = Q.shape[1] // 2
    out = np.empty(N)
    _lib.hd_batch(_p(C), _p(Q), N, n, m, 1 if both else 0, _p(out))
    return out


def hd1(c, Q, m=48, both=True):
    c = np.ascontiguousarray(c, dtype=np.float64).reshape(8)
    Q = np.ascontiguousarray(Q, dtype=np.float64).reshape(-1)
    n = len(Q) // 2
    a = ctypes.c_double(); b = ctypes.c_double(); w = ctypes.c_double()
    h = _lib.hd_one(_p(c), _p(Q), n, m, 1 if both else 0, ctypes.byref(a), ctypes.byref(b), ctypes.byref(w))
    return h, a.value, b.value, w.value


def ep(C, Q, m=32):
    C = np.ascontiguousarray(C, dtype=np.float64).reshape(-1, 8)
    N = C.shape[0]
    Q = np.ascontiguousarray(Q, dtype=np.float64).reshape(N, -1)
    n = Q.shape[1] // 2
    out = np.empty(N)
    _lib.ep_batch(_p(C), _p(Q), N, n, m, _p(out))
    return out


def resid(c, Q, m=16):
    c = np.ascontiguousarray(c, dtype=np.float64).reshape(8)
    Q = np.ascontiguousarray(Q, dtype=np.float64).reshape(-1)
    n = len(Q) // 2
    M = n * m + 1
    r = np.empty(M); gx = np.empty(M); gy = np.empty(M); tf = np.empty(M)
    _lib.resid(_p(c), _p(Q), n, m, _p(r), _p(gx), _p(gy), _p(tf))
    return r, gx, gy, tf


if __name__ == '__main__':
    import json, time
    from common import CORPUS, family_offcurves, ramp
    from geom import hausdorff_batch, cubic_at, spline_at
    C = np.array(json.load(open(CORPUS)), dtype=float)
    rng = np.random.default_rng(11)
    idx = rng.choice(len(C), 3000, replace=False)
    Cs = C[idx]

    def test_splines(Cx, n, noise, r):
        return family_offcurves(Cx, ramp(n) if n > 1 else np.array([0.0])) + r.normal(0, noise, size=(len(Cx), n, 2))

    print('1. C helper against the numpy reference geom.py (both m = 24, no refinement in numpy), 300 corpus cubics:')
    for n in (1, 2, 3, 5):
        Q = test_splines(Cs[:300], n, 0.4, rng)
        a = hd(Cs[:300], Q, m=24); b = hausdorff_batch(Cs[:300], Q, m=24)
        ok = b > 1e-3
        rel = (a - b)[ok] / b[ok]
        print(f'   n = {n}: C/numpy - 1: min {rel.min():+.1e}  median {np.median(rel):+.1e}  max {rel.max():+.1e}', flush=True)

    print('2. C helper at m samples per quadratic against itself at m = 192, 3000 corpus cubics, n = 1..5:')
    for label, noise in (('cu2qu ramp', 0.0), ('ramp + 0.4 unit noise', 0.4)):
        Qs = {n: test_splines(Cs, n, noise, np.random.default_rng(n)) for n in (1, 2, 3, 4, 5)}
        ref = {n: hd(Cs, Qs[n], m=192) for n in Qs}
        for m in (8, 16, 24, 48):
            rel = np.concatenate([(hd(Cs, Qs[n], m=m) / ref[n] - 1)[ref[n] > 1e-3] for n in Qs])
            print(f'   {label:22s} m = {m:2d}: min {rel.min():+.1e}  1% {np.quantile(rel, 0.01):+.1e}  median {np.median(rel):+.1e}'
                  f'  99% {np.quantile(rel, 0.99):+.1e}  max {rel.max():+.1e}', flush=True)

    def seg_dist(P, A, B):
        out = np.empty(len(P))
        AB = B - A
        L2 = (AB * AB).sum(-1); L2 = np.where(L2 > 0, L2, 1.0)
        for i0 in range(0, len(P), 250):
            p = P[i0:i0 + 250, None, :]
            w = np.clip(((p - A[None]) * AB[None]).sum(-1) / L2[None], 0, 1)
            proj = A[None] + w[..., None] * AB[None]
            out[i0:i0 + 250] = np.sqrt(((p - proj) ** 2).sum(-1)).min(-1)
        return out

    print('3. C helper (m = 24 and m = 192) against a polyline brute force (8001 points per curve), 40 corpus cubics:')
    rels = {24: [], 192: []}
    r2 = np.random.default_rng(5)
    for i in idx[:40]:
        c = C[i]
        n = int(r2.integers(1, 6))
        Q = test_splines(c[None], n, 0.3, r2)[0]
        ts = np.linspace(0, 1, 8001)
        Bp = cubic_at(c[None], ts[None])[0][0]
        Sp = spline_at(c[None], Q[None], (ts * n)[None])[0][0]
        hb = max(seg_dist(Sp[::2], Bp[:-1], Bp[1:]).max(), seg_dist(Bp[::2], Sp[:-1], Sp[1:]).max())
        if hb > 1e-3:
            for m in rels:
                rels[m].append(hd1(c, Q, m=m)[0] / hb - 1)
    for m in rels:
        r = np.array(rels[m])
        print(f'   m = {m:3d}: C/brute - 1: min {r.min():+.1e}  median {np.median(r):+.1e}  max {r.max():+.1e}  ({len(r)} cases)')

    Q = family_offcurves(C, ramp(4))
    t0 = time.time(); a = hd(C, Q, m=24); t1 = time.time()
    print(f'timing: {len(C)} cubics, n = 4, m = 24: {t1 - t0:.2f} s')
