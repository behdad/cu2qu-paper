#!/usr/bin/env python3
"""Shared pieces for the 'weights' research.

Notation (as in the brief).  A cubic c0 c1 c2 c3, B(t) = c0 + t b1 + t^2 b2 + t^3 b3,
b3 = d3 = c3 - 3 c2 + 3 c1 - c0.  Cut into n equal pieces; piece k has control
points p0..p3, third difference d = d3 / n^3, and "midpoint control"
q_M = (3 (p1 + p2) - p0 - p3) / 4.  The d3 family of splines has off-curve points

    Q_k = q_M,k + s_k d / 4

and implied on-curve points halfway between consecutive off-curve points.
Quadratic k minus cubic piece k at equal local parameter u is d * phi_k(u),

    phi_k(u) = (1-u)^2 g_k + u^2 g_{k+1} + u (1-u) (u - h_k)
    g_k = (s_{k-1} + s_k) / 8  (g_0 = g_n = 0),   h_k = (1 - s_k) / 2.

Running this file checks those identities on random cubics and checks the
ramp s_k = 1 - 2k/(n-1) against fontTools.
"""
import math, os, sys
import numpy as np

S = '/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad'
CORPUS = S + '/ft_cubics.json'
FT_LIB = S + '/ft-old/Lib'


# ---------------------------------------------------------------- phi algebra

def g_of(s):
    """g_0..g_n for weights s (length n)."""
    s = np.asarray(s, dtype=float)
    n = len(s)
    g = np.zeros(n + 1)
    g[1:n] = (s[:-1] + s[1:]) / 8.0
    return g

def phi_coeffs(s):
    """(n, 4) array: phi_k(u) = c0 + c1 u + c2 u^2 + c3 u^3."""
    s = np.asarray(s, dtype=float)
    n = len(s)
    g = g_of(s)
    h = (1.0 - s) / 2.0
    co = np.empty((n, 4))
    co[:, 0] = g[:-1]
    co[:, 1] = -(2 * g[:-1] + h)
    co[:, 2] = 1 + g[:-1] + g[1:] + h
    co[:, 3] = -1.0
    return co

def cubic_extrema(co):
    """Parameters in [0, 1] where c0 + c1 u + c2 u^2 + c3 u^3 can be extreme."""
    cand = [0.0, 1.0]
    a, b, c = 3 * co[3], 2 * co[2], co[1]
    if abs(a) < 1e-300:
        if abs(b) > 1e-300:
            cand.append(-c / b)
    else:
        disc = b * b - 4 * a * c
        if disc >= 0:
            r = math.sqrt(disc)
            # numerically stable roots
            q = -0.5 * (b + (r if b >= 0 else -r))
            if q != 0:
                cand.append(c / q)
            cand.append(q / a)
    return [u for u in cand if 0.0 <= u <= 1.0]

def poly(co, u):
    return co[0] + u * (co[1] + u * (co[2] + u * co[3]))

def phi_max(s, per_piece=False):
    """Exact max over pieces and u of |phi_k(u)| (critical points, not a grid)."""
    co = phi_coeffs(s)
    out = []
    for k in range(len(co)):
        best = (0.0, 0.0, 0.0)
        for u in cubic_extrema(co[k]):
            v = poly(co[k], u)
            if abs(v) > best[0]:
                best = (abs(v), u, v)
        out.append(best)
    if per_piece:
        return out
    return max(b[0] for b in out)

def ramp(n):
    return np.array([1 - 2 * k / (n - 1) for k in range(n)]) if n > 1 else np.array([0.0])

def alternating(n):
    return np.array([1.0 if k % 2 == 0 else -1.0 for k in range(n)])


# ---------------------------------------------------------------- geometry

def power_basis(c):
    """c: (..., 4, 2) -> b0, b1, b2, b3 with B(t) = b0 + t b1 + t^2 b2 + t^3 b3."""
    c = np.asarray(c, dtype=float)
    c0, c1, c2, c3 = c[..., 0, :], c[..., 1, :], c[..., 2, :], c[..., 3, :]
    return c0, 3 * (c1 - c0), 3 * (c0 - 2 * c1 + c2), c3 - 3 * c2 + 3 * c1 - c0

def qM_points(c, n):
    """q_M of the n equal pieces: (..., n, 2).  Uses
    q_M = (B(t0) + B(t1)) / 2 - (t1 - t0) (B'(t1) - B'(t0)) / 4."""
    b0, b1, b2, b3 = power_basis(c)
    k = np.arange(n)
    t0 = (k / n)[:, None]
    t1 = ((k + 1) / n)[:, None]
    def B(t):
        return b0[..., None, :] + t * (b1[..., None, :] + t * (b2[..., None, :] + t * b3[..., None, :]))
    def D(t):
        return b1[..., None, :] + t * (2 * b2[..., None, :] + 3 * t * b3[..., None, :])
    return (B(t0) + B(t1)) / 2 - (D(t1) - D(t0)) / (4.0 * n)

def family_offcurves(c, s):
    """Off-curve points Q_k = q_M,k + s_k d / 4 of the d3-family spline: (..., n, 2)."""
    s = np.asarray(s, dtype=float)
    n = len(s)
    b0, b1, b2, b3 = power_basis(c)
    return qM_points(c, n) + s[:, None] * b3[..., None, :] / (4.0 * n ** 3)

def oncurves(c, Q):
    """On-curve points of the spline with off-curve points Q (..., n, 2): (..., n+1, 2)."""
    c = np.asarray(c, dtype=float)
    mids = (Q[..., :-1, :] + Q[..., 1:, :]) / 2
    return np.concatenate([c[..., 0:1, :], mids, c[..., 3:4, :]], axis=-2)


if __name__ == '__main__':
    rng = np.random.default_rng(1)
    # 1. E_k(u) = d phi_k(u) for random weights on random cubics.
    worst = 0.0
    for trial in range(3000):
        c = rng.uniform(-1000, 1000, size=(4, 2))
        n = int(rng.integers(1, 25))
        s = rng.uniform(-2, 2, size=n)
        Q = family_offcurves(c, s)
        P = oncurves(c, Q)
        b0, b1, b2, b3 = power_basis(c)
        co = phi_coeffs(s)
        d = b3 / n ** 3
        for k in range(n):
            u = np.linspace(0, 1, 9)[:, None]
            t = (k + u) / n
            cub = b0 + t * (b1 + t * (b2 + t * b3))
            quad = (1 - u) ** 2 * P[k] + 2 * u * (1 - u) * Q[k] + u * u * P[k + 1]
            pred = d * poly(co[k], u)
            worst = max(worst, np.abs(quad - cub - pred).max() / max(1.0, np.abs(b3).max()))
    print(f'identity S - B = d phi_k(u), random weights, n = 1..24: worst relative discrepancy {worst:.2e}')

    # 2. s = +1 / -1 are the points on the piece's start / end tangent.
    worst = 0.0
    for trial in range(1000):
        c = rng.uniform(-1000, 1000, size=(4, 2))
        n = int(rng.integers(1, 12))
        b0, b1, b2, b3 = power_basis(c)
        Bf = lambda t: b0 + t * (b1 + t * (b2 + t * b3))
        Df = lambda t: b1 + t * (2 * b2 + 3 * t * b3)
        Qp = family_offcurves(c, np.ones(n)); Qm = family_offcurves(c, -np.ones(n))
        for k in range(n):
            t0, t1 = k / n, (k + 1) / n
            worst = max(worst, np.abs(Qp[k] - (Bf(t0) + Df(t0) / (2 * n))).max(),
                        np.abs(Qm[k] - (Bf(t1) - Df(t1) / (2 * n))).max())
    print(f's = +1 is B(t0) + B\'(t0) l/2 and s = -1 is B(t1) - B\'(t1) l/2: worst gap {worst:.2e}')

    # 3. the ramp is fontTools' spline.
    sys.path.insert(0, FT_LIB)
    from fontTools.cu2qu.cu2qu import cubic_approx_spline
    worst = 0.0; count = 0
    for trial in range(20000):
        c = rng.uniform(-1000, 1000, size=(4, 2))
        n = int(rng.integers(2, 20))
        cc = [complex(*p) for p in c]
        sp = cubic_approx_spline(cc, n, 1e30, True)
        Q = family_offcurves(c, ramp(n))
        ft = np.array([[z.real, z.imag] for z in sp[1:-1]])
        worst = max(worst, np.abs(ft - Q).max()); count += 1
    print(f'ramp vs fontTools cubic_approx_spline on {count} random cubics (n = 2..19): worst gap {worst:.2e}')
