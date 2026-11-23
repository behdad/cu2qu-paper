#!/usr/bin/env python3
"""Tangent-intersection ("T") splines: closed-form geometry and a closed-form
first-order estimate of their geometric error.

Piece k of the cubic, between knots t_k and t_{k+1}: end points A, B, unit
tangents there, tangent intersection Q.  With a = A - Q and b = B - Q the
cubic piece is, in the oblique frame (Q; a, b),

    x(u) = (1-u)^2 (1 + ea u),      ea = 2 - l |B'(t_k)|   / |a|
    y(u) = u^2 (1 + eb (1-u)),      eb = 2 - l |B'(t_k+1)| / |b|

(l the parameter length of the piece) and a quadratic with control point Q that
touches the two tangent lines at Q + xa a and Q + xb b is  sqrt(x/xa) + sqrt(y/xb) = 1.
So the cubic's point at u is off that parabola by, to first order,

    D(u) = F / |grad F|,   F = sqrt(x/xa) + sqrt(y/xb) - 1,
    |grad F| = | b / (2 sqrt(x xa)) - a / (2 sqrt(y xb)) | / |a x b|.

For the T spline with implied on-curve points, xa and xb say where the
midpoints of consecutive Q's fall on the tangent lines (1 = on the cubic).
For one quadratic (xa = xb = 1), small ea + eb:  max D = (ea + eb) h / 8, h the
height of the triangle A Q B over its base.

    t_pieces(C, knots)            geometry of every piece
    t_estimate(C, knots)          estimated Hausdorff distance of the T spline
"""
import numpy as np
from common import power_basis
from schemes import unit_end_tangents, cross


def t_pieces(C, knots):
    """Per piece (N, n): dict with A, B, Q, a, b, ea, eb, valid, and the spline's off-curve points."""
    C = np.asarray(C, dtype=float)
    N = len(C); n = knots.shape[1] - 1
    b0, b1, b2, b3 = power_basis(C)
    t = knots[:, :, None]
    P = b0[:, None] + t * (b1[:, None] + t * (b2[:, None] + t * b3[:, None]))
    Dv = b1[:, None] + t * (2 * b2[:, None] + 3 * t * b3[:, None])
    speed = np.hypot(Dv[..., 0], Dv[..., 1])
    e0, e1 = unit_end_tangents(C)
    Dir = Dv / np.where(speed > 1e-12, speed, 1.0)[..., None]
    Dir[:, 0] = e0; Dir[:, -1] = e1
    A = P[:, :-1]; B = P[:, 1:]; Da = Dir[:, :-1]; Db = Dir[:, 1:]
    den = cross(Da, Db)
    chord = B - A
    ok = np.abs(den) > 1e-9
    den_s = np.where(ok, den, 1.0)
    la = cross(chord, Db) / den_s            # Q = A + la Da
    lb = cross(Da, chord) / den_s            # Q = B - lb Db
    cl = np.hypot(chord[..., 0], chord[..., 1])
    valid = ok & (la > 1e-9) & (lb > 1e-9) & (la < 4 * cl) & (lb < 4 * cl)
    la_s = np.where(valid, la, 1.0); lb_s = np.where(valid, lb, 1.0)
    Q = A + la_s[..., None] * Da
    l = knots[:, 1:] - knots[:, :-1]
    ea = 2 - l * speed[:, :-1] / la_s
    eb = 2 - l * speed[:, 1:] / lb_s
    valid = valid & (ea > -0.999) & (eb > -0.999)
    return dict(A=A, B=B, Q=Q, a=A - Q, b=B - Q, ea=ea, eb=eb, valid=valid, la=la_s, lb=lb_s, l=l)


def t_profile(pc, xa, xb, U=33):
    """D(u) on a grid, shape (N, n, U); xa, xb (N, n)."""
    u = np.linspace(0.0, 1.0, U)[None, None, :]
    ea = pc['ea'][..., None]; eb = pc['eb'][..., None]
    x = (1 - u) ** 2 * (1 + ea * u)
    y = u ** 2 * (1 + eb * (1 - u))
    xa = xa[..., None]; xb = xb[..., None]
    sx = np.sqrt(np.maximum(x * xa, 0.0)); sy = np.sqrt(np.maximum(y * xb, 0.0))
    F = sx / xa + sy / xb - 1
    a = pc['a'][:, :, None, :]; b = pc['b'][:, :, None, :]
    axb = np.abs(cross(pc['a'], pc['b']))[..., None]
    v = b * sy[..., None] - a * sx[..., None]
    den = np.hypot(v[..., 0], v[..., 1])
    return F * 2 * axb * sx * sy / np.where(den > 1e-300, den, 1.0)


def t_estimate(C, knots, U=33, detail=False):
    """Estimated geometric error of the T spline at the given knots: (N,), NaN where a piece has no valid T point."""
    pc = t_pieces(C, knots)
    Q = pc['Q']; n = Q.shape[1]
    # where the implied on-curve points fall on the tangent lines
    xa = np.ones(Q.shape[:2]); xb = np.ones(Q.shape[:2])
    if n > 1:
        M = (Q[:, :-1] + Q[:, 1:]) / 2                       # on-curve point between pieces k and k+1
        a = pc['a']; b = pc['b']
        # piece k+1 starts at M[k]:  M = Q_{k+1} + xa a_{k+1}
        xa[:, 1:] = ((M - Q[:, 1:]) * a[:, 1:]).sum(-1) / (a[:, 1:] ** 2).sum(-1)
        # piece k ends at M[k]:      M = Q_k + xb b_k
        xb[:, :-1] = ((M - Q[:, :-1]) * b[:, :-1]).sum(-1) / (b[:, :-1] ** 2).sum(-1)
    ok = pc['valid'].all(1) & (xa > 1e-6).all(1) & (xb > 1e-6).all(1)
    D = t_profile(pc, np.where(xa > 1e-6, xa, 1.0), np.where(xb > 1e-6, xb, 1.0), U)
    est = np.abs(D).max(axis=(1, 2))
    est = np.where(ok, est, np.nan)
    if detail:
        return est, D, pc, xa, xb
    return est


if __name__ == '__main__':
    import json, os
    import schemes
    from common import CORPUS
    BIG = '/home/behdad/hb/cu2qu-work/research-weights'
    C = np.array(json.load(open(CORPUS)), dtype=float)
    print('closed-form estimate against the measured Hausdorff distance, T spline, whole corpus')
    print('(cubics whose pieces all have a valid tangent intersection and whose measured distance is > 0.01):')
    print('   knots   n   cubics   estimate/measured: 1%      10%     median   90%     99%     max   | share with all pieces valid')
    for label, key, knots_of in (('equal t', 'Tuni', lambda n: schemes.uniform_knots(len(C), n)),
                                 ('gamma  ', 'Tgam', lambda n: schemes.gamma_knots(C, n))):
        H = np.load(os.path.join(BIG, f'corpus_{key}.npy'))
        for n in (1, 2, 3, 4):
            est = t_estimate(C, knots_of(n))
            h = H[:, n - 1]
            ok = np.isfinite(est) & np.isfinite(h) & (h > 0.01)
            r = est[ok] / h[ok]
            q = np.quantile(r, [0.01, 0.1, 0.5, 0.9, 0.99])
            print(f'   {label} {n:2d}  {ok.sum():6d}   ' + '  '.join(f'{v:6.3f}' for v in q) + f'  {r.max():6.3f}  |  {np.isfinite(est).mean():.3f}')
