#!/usr/bin/env python3
"""Geometric distance between a cubic and an implied-on-curve quadratic spline.

A measuring instrument, not a conversion algorithm: dense sampling plus a few
Newton steps for the nearest point.

    hausdorff(c, Q)          one cubic c (4, 2), off-curve points Q (n, 2)
    hausdorff_batch(C, Q)    C (N, 4, 2), Q (N, n, 2), chunked
    equal_param_batch(C, Q)  the cu2qu measure: largest |S((k+u)) - B((k+u)/n)|

Both directed distances are computed (spline -> cubic and cubic -> spline);
the Hausdorff distance is the larger.  The spline is S(xi), xi in [0, n],
quadratic k on [k, k+1]: on_k (1-u)^2 + 2 Q_k u (1-u) + on_{k+1} u^2 with
on_0 = c0, on_n = c3 and on_k = (Q_{k-1} + Q_k) / 2 otherwise.
"""
import numpy as np


def _power(C):
    c0, c1, c2, c3 = C[..., 0, :], C[..., 1, :], C[..., 2, :], C[..., 3, :]
    return c0, 3 * (c1 - c0), 3 * (c0 - 2 * c1 + c2), c3 - 3 * c2 + 3 * c1 - c0


def cubic_at(C, t):
    """C (N, 4, 2), t (N, M) -> points, first and second derivatives (N, M, 2)."""
    b0, b1, b2, b3 = _power(C)
    t = t[..., None]
    b0 = b0[:, None, :]; b1 = b1[:, None, :]; b2 = b2[:, None, :]; b3 = b3[:, None, :]
    P = b0 + t * (b1 + t * (b2 + t * b3))
    D = b1 + t * (2 * b2 + 3 * t * b3)
    DD = 2 * b2 + 6 * t * b3
    return P, D, DD


def spline_ctrl(C, Q):
    """On-curve points (N, n+1, 2) of the spline."""
    mids = (Q[:, :-1] + Q[:, 1:]) / 2
    return np.concatenate([C[:, 0:1], mids, C[:, 3:4]], axis=1)


def spline_at(C, Q, xi):
    """xi (N, M) in [0, n] -> points, first and second derivatives (N, M, 2)."""
    n = Q.shape[1]
    on = spline_ctrl(C, Q)
    k = np.clip(np.floor(xi).astype(int), 0, n - 1)
    u = (xi - k)[..., None]
    idx = np.arange(Q.shape[0])[:, None]
    a = on[idx, k]; q = Q[idx, k]; b = on[idx, k + 1]
    P = (1 - u) ** 2 * a + 2 * u * (1 - u) * q + u * u * b
    D = 2 * (1 - u) * (q - a) + 2 * u * (b - q)
    DD = 2 * (a - 2 * q + b)
    return P, D, DD


def _nearest_on_cubic(C, P, K, iters=5):
    """For points P (N, M, 2): parameter of the nearest point of the cubic, and distance."""
    N, M = P.shape[:2]
    tj = np.linspace(0.0, 1.0, K)
    Bj, _, _ = cubic_at(C, np.broadcast_to(tj, (N, K)))
    d2 = ((P[:, :, None, :] - Bj[:, None, :, :]) ** 2).sum(-1)       # (N, M, K)
    j = d2.argmin(-1)
    t0 = tj[j]
    best = np.sqrt(np.take_along_axis(d2, j[..., None], -1)[..., 0])
    lo = np.clip(t0 - 1.0 / (K - 1), 0, 1); hi = np.clip(t0 + 1.0 / (K - 1), 0, 1)
    t = t0.copy()
    for _ in range(iters):
        B, D, DD = cubic_at(C, t)
        v = B - P
        f = (v * D).sum(-1)
        fp = (D * D).sum(-1) + (v * DD).sum(-1)
        step = np.where(fp > 1e-30, f / np.where(fp > 1e-30, fp, 1.0), np.sign(f) * (hi - lo) * 0.25)
        t = np.clip(t - step, lo, hi)
    B, D, _ = cubic_at(C, t)
    dist = np.sqrt(((B - P) ** 2).sum(-1))
    worse = dist > best
    t = np.where(worse, t0, t); dist = np.where(worse, best, dist)
    return t, dist


def _nearest_on_spline(C, Q, P, m, iters=5):
    """For points P (N, M, 2): parameter xi of the nearest point of the spline, and distance."""
    N, n = Q.shape[:2]
    K = n * m + 1
    xj = np.linspace(0.0, n, K)
    Sj, _, _ = spline_at(C, Q, np.broadcast_to(xj, (N, K)))
    d2 = ((P[:, :, None, :] - Sj[:, None, :, :]) ** 2).sum(-1)
    j = d2.argmin(-1)
    x0 = xj[j]
    best = np.sqrt(np.take_along_axis(d2, j[..., None], -1)[..., 0])
    lo = np.clip(x0 - 1.0 / m, 0, n); hi = np.clip(x0 + 1.0 / m, 0, n)
    x = x0.copy()
    for _ in range(iters):
        Sx, D, DD = spline_at(C, Q, x)
        v = Sx - P
        f = (v * D).sum(-1)
        fp = (D * D).sum(-1) + (v * DD).sum(-1)
        step = np.where(fp > 1e-30, f / np.where(fp > 1e-30, fp, 1.0), np.sign(f) * (hi - lo) * 0.25)
        x = np.clip(x - step, lo, hi)
    Sx, _, _ = spline_at(C, Q, x)
    dist = np.sqrt(((Sx - P) ** 2).sum(-1))
    worse = dist > best
    x = np.where(worse, x0, x); dist = np.where(worse, best, dist)
    return x, dist


def hausdorff_chunk(C, Q, m=24, both=True, detail=False):
    """C (N, 4, 2), Q (N, n, 2) -> (N,) Hausdorff distance (or the two directed ones)."""
    N, n = Q.shape[:2]
    M = n * m + 1
    xi = np.broadcast_to(np.linspace(0.0, n, M), (N, M))
    Sp, _, _ = spline_at(C, Q, xi)
    K = max(4 * n * m // 3, 48) + 1
    t1, d1 = _nearest_on_cubic(C, Sp, K)
    a = d1.max(-1)
    if not both:
        return a
    tt = np.broadcast_to(np.linspace(0.0, 1.0, M), (N, M))
    Bp, _, _ = cubic_at(C, tt)
    x2, d2 = _nearest_on_spline(C, Q, Bp, m)
    b = d2.max(-1)
    if detail:
        return a, b, (xi, d1, t1), (tt, d2, x2)
    return np.maximum(a, b)


def hausdorff_batch(C, Q, m=24, both=True, chunk=None):
    C = np.asarray(C, dtype=float); Q = np.asarray(Q, dtype=float)
    N, n = Q.shape[:2]
    if chunk is None:
        chunk = max(1, int(4e6 / ((n * m + 1) * (max(4 * n * m // 3, 48) + 1))))
    out = np.empty(N)
    for i in range(0, N, chunk):
        out[i:i + chunk] = hausdorff_chunk(C[i:i + chunk], Q[i:i + chunk], m, both)
    return out


def hausdorff(c, Q, m=48, both=True):
    return float(hausdorff_chunk(np.asarray(c, dtype=float)[None], np.asarray(Q, dtype=float)[None], m, both)[0])


def equal_param_batch(C, Q, m=32):
    """Largest distance between spline and cubic at equal parameters (uniform pieces)."""
    C = np.asarray(C, dtype=float); Q = np.asarray(Q, dtype=float)
    N, n = Q.shape[:2]
    M = n * m + 1
    xi = np.broadcast_to(np.linspace(0.0, n, M), (N, M))
    Sp, _, _ = spline_at(C, Q, xi)
    Bp, _, _ = cubic_at(C, xi / n)
    return np.sqrt(((Sp - Bp) ** 2).sum(-1)).max(-1)


if __name__ == '__main__':
    # self-test against brute force on random cubics and perturbed splines
    from common import family_offcurves, ramp
    rng = np.random.default_rng(3)
    worst = 0.0
    for trial in range(300):
        c = rng.uniform(-500, 500, size=(4, 2))
        n = int(rng.integers(2, 7))
        Q = family_offcurves(c, ramp(n)) + rng.normal(0, 3.0, size=(n, 2))
        h = hausdorff(c, Q, m=48)
        # brute force: 4000 samples on each curve
        ts = np.linspace(0, 1, 4001)
        Bp, _, _ = cubic_at(c[None], ts[None]); Bp = Bp[0]
        Sp, _, _ = spline_at(c[None], Q[None], (ts * n)[None]); Sp = Sp[0]
        D = np.sqrt(((Bp[:, None, :] - Sp[None, :, :]) ** 2).sum(-1))
        hb = max(D.min(0).max(), D.min(1).max())
        worst = max(worst, abs(h - hb) / max(hb, 1e-9))
    print(f'hausdorff() vs 4001x4001 brute force on 300 random cubics: worst relative difference {worst:.2e}')
