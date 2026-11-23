#!/usr/bin/env python3
"""Closed-form spline schemes, vectorised over cubics C (N, 4, 2).

Each returns off-curve points Q (N, n, 2); the on-curve points are implied.

d3 family on n equal pieces (Q_k = q_M,k + s_k d / 4):
  ramp(C, n)        cu2qu:  s_k = 1 - 2k/(n-1)
  optu(C, n)        least equal-parameter error with both end tangents kept:
                    alternating +1, -1 for even n; the LP optimum for odd n
  zero(C, n)        s = 0: least equal-parameter error, end tangents not kept
  relaxed(C, n, x)  +x, -x, ... for even n; LP optimum with s_0 = x for odd n

Outside the d3 family:
  tangents(C, knots)  off-curve point k is where the cubic's tangents at knots
                      k and k+1 meet ("T"); keeps both end tangents; a piece
                      whose tangents do not meet ahead of both ends falls back
                      to the d3-family point (s = +1 at the start piece, -1 at
                      the end piece, 0 inside)
  uniform_knots(N, n), affine_knots(C, n), gamma_knots(C, n)
"""
import math
import numpy as np
from common import family_offcurves, ramp as ramp_w, alternating, power_basis

_OPT = {}


def optu_weights(n, x=1.0):
    """Equal-parameter-optimal weights with s_0 = x, s_{n-1} = -x (antisymmetric)."""
    key = (n, x)
    if key not in _OPT:
        if n == 1:
            _OPT[key] = np.array([0.0])
        elif n % 2 == 0:
            _OPT[key] = x * alternating(n)
        else:
            from q1_lp import solve
            t, s = solve(n, f'relax:{x}', antisym=True, secondary='l1')
            _OPT[key] = s
    return _OPT[key]


def ramp(C, n):
    return family_offcurves(C, ramp_w(n) if n > 1 else np.array([0.0]))

def optu(C, n):
    return family_offcurves(C, optu_weights(n))

def zero(C, n):
    return family_offcurves(C, np.zeros(n))

def relaxed(C, n, x):
    return family_offcurves(C, optu_weights(n, x))


def cross(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def unit_end_tangents(C):
    """(N, 2) unit forward tangents at t = 0 and t = 1, with cu2qu-like fallbacks."""
    C = np.asarray(C, dtype=float)
    def pick(cands):
        out = np.zeros_like(cands[0]); done = np.zeros(len(out), dtype=bool)
        for v in cands:
            l = np.hypot(v[:, 0], v[:, 1])
            use = (~done) & (l > 1e-12)
            out[use] = v[use] / l[use][:, None]
            done |= use
        out[~done] = np.array([1.0, 0.0])
        return out
    e0 = pick([C[:, 1] - C[:, 0], C[:, 2] - C[:, 0], C[:, 3] - C[:, 0]])
    e1 = pick([C[:, 3] - C[:, 2], C[:, 3] - C[:, 1], C[:, 3] - C[:, 0]])
    return e0, e1


def uniform_knots(N, n):
    return np.broadcast_to(np.linspace(0.0, 1.0, n + 1), (N, n + 1)).copy()


def curvature_poly(C):
    """A, Cc, D with  B' x B'' = 18 (A (1-t)^2 + Cc t (1-t) + D t^2)."""
    C = np.asarray(C, dtype=float)
    e0 = C[:, 1] - C[:, 0]; e1 = C[:, 2] - C[:, 1]; e2 = C[:, 3] - C[:, 2]
    return cross(e0, e1), cross(e0, e2), cross(e1, e2)


def affine_knots(C, n, grid=1025):
    """Knots at equal affine arc length, int |B' x B''|^(1/3) dt (numerical)."""
    A, Cc, D = curvature_poly(C)
    t = np.linspace(0.0, 1.0, grid)
    p = A[:, None] * (1 - t) ** 2 + Cc[:, None] * t * (1 - t) + D[:, None] * t * t
    w = np.cbrt(np.abs(p))
    cum = np.concatenate([np.zeros((len(w), 1)), np.cumsum((w[:, 1:] + w[:, :-1]) / 2, axis=1) / (grid - 1)], axis=1)
    tot = cum[:, -1]
    out = uniform_knots(len(w), n)
    for i in range(len(w)):
        if tot[i] > 1e-9:
            out[i] = np.interp(np.linspace(0, tot[i], n + 1), cum[i], t)
    out[:, 0] = 0.0; out[:, -1] = 1.0
    return out


def gamma_knots(C, n, clamp=0.6):
    """Closed-form stand-in for affine_knots: t = x + g x (1 - x) at x = k/n, with
    (1 + g) / (1 - g) = (|D| / |A|)^(1/3): the knot density is right at both ends."""
    A, Cc, D = curvature_poly(C)
    a = np.cbrt(np.abs(A)); d = np.cbrt(np.abs(D))
    g = np.where(a + d > 1e-12, (d - a) / np.where(a + d > 1e-12, a + d, 1.0), 0.0)
    g = np.clip(g, -clamp, clamp)
    x = np.linspace(0.0, 1.0, n + 1)[None, :]
    return x + g[:, None] * x * (1 - x)


def tangents(C, knots, return_valid=False):
    """Off-curve points at the intersections of consecutive knot tangents."""
    C = np.asarray(C, dtype=float)
    N = len(C); n = knots.shape[1] - 1
    b0, b1, b2, b3 = power_basis(C)
    t = knots[:, :, None]
    P = b0[:, None] + t * (b1[:, None] + t * (b2[:, None] + t * b3[:, None]))
    Dv = b1[:, None] + t * (2 * b2[:, None] + 3 * t * b3[:, None])
    e0, e1 = unit_end_tangents(C)
    # directions: unit vectors (the ends may have a zero derivative)
    L = np.hypot(Dv[..., 0], Dv[..., 1])
    Dir = Dv / np.where(L > 1e-12, L, 1.0)[..., None]
    Dir[:, 0] = e0; Dir[:, -1] = e1
    Pa = P[:, :-1]; Pb = P[:, 1:]; Da = Dir[:, :-1]; Db = Dir[:, 1:]
    den = cross(Da, Db)
    chord = Pb - Pa
    ok = np.abs(den) > 1e-9
    den_s = np.where(ok, den, 1.0)
    a = cross(chord, Db) / den_s            # Pa + a Da
    b = cross(Da, chord) / den_s            # = Pb - b Db
    cl = np.hypot(chord[..., 0], chord[..., 1])
    valid = ok & (a > 1e-9) & (b > 1e-9) & (a < 4 * cl) & (b < 4 * cl)
    Q = Pa + a[..., None] * Da
    # fallback: the d3-family point of the (possibly unequal) piece
    l = (knots[:, 1:] - knots[:, :-1])[..., None]
    Dva = Dv[:, :-1]; Dvb = Dv[:, 1:]
    qm = (Pa + Pb) / 2 - l * (Dvb - Dva) / 4
    d = b3[:, None] * l ** 3
    fb = qm.copy()
    fb[:, 0] = qm[:, 0] + d[:, 0] / 4          # on the start tangent
    if n > 1:
        fb[:, -1] = qm[:, -1] - d[:, -1] / 4   # on the end tangent
    Q = np.where(valid[..., None], Q, fb)
    if return_valid:
        return Q, valid
    return Q


def inflections(C):
    """Number of sign changes of B' x B'' = 18 (A (1-t)^2 + Cc t (1-t) + D t^2) on (0, 1)."""
    A, Cc, D = curvature_poly(C)
    qa = A - Cc + D; qb = Cc - 2 * A; qc = A
    cnt = np.zeros(len(A), dtype=int)
    scale = np.abs(A) + np.abs(Cc) + np.abs(D) + 1e-300
    disc = qb * qb - 4 * qa * qc
    lin = np.abs(qa) < 1e-12 * scale
    with np.errstate(all='ignore'):
        r_lin = -qc / qb
        cnt += (lin & (np.abs(qb) > 1e-12 * scale) & (r_lin > 1e-9) & (r_lin < 1 - 1e-9)).astype(int)
        sq = np.sqrt(np.maximum(disc, 0))
        for sgn in (1, -1):
            r = (-qb + sgn * sq) / (2 * qa)
            cnt += (~lin & (disc > 1e-18 * scale ** 2) & (r > 1e-9) & (r < 1 - 1e-9)).astype(int)
    return cnt
