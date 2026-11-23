#!/usr/bin/env python3
"""Numerical minimax fit of an implied-on-curve quadratic spline to a cubic,
under the geometric distance.  A yardstick, not a conversion algorithm.

    best_spline(c, n, mode, starts)  ->  (hausdorff, Q)

Families (mode), all linear maps x -> off-curve points Q = base + L x:
  'tan'    both end tangents kept: Q_0 on the ray c0 + a e0, Q_{n-1} on the ray
           c3 - b e1, interior off-curve points free          (2n - 2 unknowns)
  'free'   every off-curve point free                         (2n unknowns)
  'w'      d3 family on equal pieces, Q_k = q_M,k + s_k d/4,
           s_0 = +1, s_{n-1} = -1                             (n - 2 unknowns)
  'ws'     as 'w', and the two end points may slide along the
           end tangents                                       (n unknowns)
  'w0'     d3 family, all weights free (tangents not kept)    (n unknowns)

Method: sequential linear programming on the signed distances r_i from spline
samples to the cubic (from chd.resid): minimise t subject to
|r_i + grad_i . dQ_i| <= t inside a trust region, repeat.  The value reported
is the two-sided Hausdorff distance from chd.hd1, not the LP's.
"""
import math
import numpy as np
from scipy.optimize import linprog
import chd
from common import power_basis, qM_points

HIGHS = {'primal_feasibility_tolerance': 1e-9, 'dual_feasibility_tolerance': 1e-9}


def end_tangents(c):
    """Unit vectors along the cubic at its start and its end (both pointing forward)."""
    c = np.asarray(c, dtype=float)
    def unit(v):
        l = math.hypot(v[0], v[1])
        return v / l if l > 1e-12 else None
    e0 = unit(c[1] - c[0])
    if e0 is None: e0 = unit(c[2] - c[0])
    if e0 is None: e0 = unit(c[3] - c[0])
    e1 = unit(c[3] - c[2])
    if e1 is None: e1 = unit(c[3] - c[1])
    if e1 is None: e1 = unit(c[3] - c[0])
    if e0 is None: e0 = np.array([1.0, 0.0])
    if e1 is None: e1 = np.array([1.0, 0.0])
    return e0, e1


def weight_matrix(n, m):
    """Spline point at xi_i = i/m is  W[i] . Q + w0[i] c0 + w3[i] c3."""
    M = n * m + 1
    W = np.zeros((M, n)); w0 = np.zeros(M); w3 = np.zeros(M)
    for i in range(M):
        k = min(i // m, n - 1)
        u = i / m - k
        W[i, k] += 2 * u * (1 - u)
        if k >= 1:
            W[i, k - 1] += (1 - u) ** 2 / 2; W[i, k] += (1 - u) ** 2 / 2
        else:
            w0[i] += (1 - u) ** 2
        if k <= n - 2:
            W[i, k] += u * u / 2; W[i, k + 1] += u * u / 2
        else:
            w3[i] += u * u
    return W, w0, w3


_WM = {}


def family(c, n, mode):
    """(base (n, 2), L (n, 2, p), lower bounds (p,) or None)."""
    c = np.asarray(c, dtype=float)
    e0, e1 = end_tangents(c)
    b0, b1, b2, b3 = power_basis(c)
    d = b3 / n ** 3
    dl = math.hypot(d[0], d[1])
    dh = d / dl if dl > 1e-12 else np.array([1.0, 0.0])
    qm = qM_points(c, n)
    if mode == 'free':
        p = 2 * n
        base = np.zeros((n, 2)); L = np.zeros((n, 2, p))
        for k in range(n):
            L[k, 0, 2 * k] = 1; L[k, 1, 2 * k + 1] = 1
        return base, L, None
    if mode == 'tan':
        if n == 1:
            return None
        p = 2 * (n - 2) + 2
        base = np.zeros((n, 2)); L = np.zeros((n, 2, p))
        base[0] = c[0]; L[0, :, 0] = e0
        base[n - 1] = c[3]; L[n - 1, :, p - 1] = -e1
        for k in range(1, n - 1):
            L[k, 0, 1 + 2 * (k - 1)] = 1; L[k, 1, 2 + 2 * (k - 1)] = 1
        lb = np.full(p, -np.inf); lb[0] = 1e-6; lb[p - 1] = 1e-6
        return base, L, lb
    if mode in ('w', 'ws', 'w0'):
        # unknowns in font units: the displacement of Q_k along the unit vector of d
        if mode == 'w':
            if n < 3:
                return qm + np.array([1.0, -1.0][:n])[:, None] * d / 4 if n == 2 else None, np.zeros((n, 2, 0)), None
            p = n - 2
            base = qm.copy(); L = np.zeros((n, 2, p))
            base[0] += d / 4; base[n - 1] -= d / 4
            for k in range(1, n - 1):
                L[k, :, k - 1] = dh
            return base, L, None
        if mode == 'ws':
            if n < 2:
                return None
            p = n
            base = qm.copy(); L = np.zeros((n, 2, p))
            base[0] += d / 4; base[n - 1] -= d / 4
            L[0, :, 0] = e0; L[n - 1, :, n - 1] = -e1
            for k in range(1, n - 1):
                L[k, :, k] = dh
            return base, L, None
        p = n
        base = qm.copy(); L = np.zeros((n, 2, p))
        for k in range(n):
            L[k, :, k] = dh
        return base, L, None
    raise ValueError(mode)


def project(base, L, Q):
    """Least-squares x with base + L x ~ Q."""
    n, _, p = L.shape
    if p == 0:
        return np.zeros(0)
    A = L.reshape(2 * n, p)
    x, *_ = np.linalg.lstsq(A, (Q - base).reshape(2 * n), rcond=None)
    return x


def slp(c, n, base, L, lb, x0, m=16, iters=60, verbose=False):
    """Minimise max_i |r_i| over x.  Returns (x, value at the m-sampling)."""
    c = np.asarray(c, dtype=float)
    p = L.shape[2]
    key = (n, m)
    if key not in _WM:
        _WM[key] = weight_matrix(n, m)
    W, w0, w3 = _WM[key]
    Lx = W @ L[:, 0, :]          # (M, p): d(P_i.x)/dx
    Ly = W @ L[:, 1, :]

    def F(x):
        Q = base + (L @ x)
        r, gx, gy, tf = chd.resid(c, Q, m)
        return np.abs(r).max(), r, gx, gy

    x = np.array(x0, dtype=float)
    if lb is not None:
        x = np.maximum(x, lb)
    f, r, gx, gy = F(x)
    if p == 0:
        return x, f
    scale = max(np.abs(c[3] - c[0]).max(), np.abs(c[1] - c[0]).max(), np.abs(c[2] - c[3]).max(), 1.0)
    delta = max(min(4 * f, 0.05 * scale), 1e-3)
    cost = np.zeros(p + 1); cost[-1] = 1.0
    for it in range(iters):
        A = gx[:, None] * Lx + gy[:, None] * Ly             # (M, p)
        Aub = np.vstack([np.hstack([A, -np.ones((len(r), 1))]),
                         np.hstack([-A, -np.ones((len(r), 1))])])
        bub = np.concatenate([-r, r])
        lo = -np.full(p, delta); hi = np.full(p, delta)
        if lb is not None:
            lo = np.maximum(lo, lb - x)
        bounds = [(lo[j], hi[j]) for j in range(p)] + [(0, None)]
        res = linprog(cost, A_ub=Aub, b_ub=bub, bounds=bounds, method='highs', options=HIGHS)
        if not res.success:
            delta *= 0.25
            if delta < 1e-9: break
            continue
        dx = res.x[:p]; pred = res.x[-1]
        f2, r2, gx2, gy2 = F(x + dx)
        if verbose:
            print(f'    it {it}: f {f:.6f} pred {pred:.6f} got {f2:.6f} delta {delta:.4g}')
        if f2 < f:
            good = (f - f2) > 0.5 * (f - pred)
            gain = f - f2
            x = x + dx; f, r, gx, gy = f2, r2, gx2, gy2
            if good and np.abs(dx).max() > 0.9 * delta:
                delta *= 2
            if gain < 1e-7 * max(f, 1e-9) and f - pred < 1e-6 * max(f, 1e-9):
                break
        else:
            delta *= 0.3
            if delta < 1e-8 or (f - pred) < 1e-7 * max(f, 1e-9):
                break
    return x, f


def best_spline(c, n, mode, starts, m=None, m_final=48, verbose=False):
    """starts: list of (n, 2) off-curve arrays to start from.  Returns (hausdorff, Q, label index)."""
    fam = family(c, n, mode)
    if fam is None:
        return None
    base, L, lb = fam
    if m is None:
        m = 24 if n <= 2 else 16
    best = None
    for si, Q0 in enumerate(starts):
        if Q0 is None or not np.all(np.isfinite(Q0)):
            continue
        x0 = project(base, L, np.asarray(Q0, dtype=float))
        if lb is not None:
            x0 = np.maximum(x0, lb)
        # the start itself (projected into the family) is a candidate too: the
        # optimiser works on the one-sided distance at a coarser sampling
        Qs = base + (L @ x0)
        h = chd.hd1(c, Qs, m=m_final)[0]
        if best is None or h < best[0]:
            best = (h, Qs, si)
        x, f = slp(c, n, base, L, lb, x0, m=m, verbose=verbose)
        Q = base + (L @ x)
        h = chd.hd1(c, Q, m=m_final)[0]
        if best is None or h < best[0]:
            best = (h, Q, si)
    return best
