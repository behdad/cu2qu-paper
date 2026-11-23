#!/usr/bin/env python3
"""A closed-form first-order estimate of the geometric distance between a cubic
and ANY implied-on-curve quadratic spline, evaluated at fixed nodes.

Quadratic k has end points P0, P1 and control point Q.  In the oblique frame
(Q; a = P0 - Q, b = P1 - Q) its parabola is  G(x, y) = (1 - x - y)^2 - 4 x y = 0.
A point of the cubic with frame coordinates (x, y) is off the parabola by

    |G| / |grad G|,   |grad G| = |G_x b - G_y a| / |a x b|,
    G_x = 2x - 2y - 2,  G_y = 2y - 2x - 2

to first order.  The cubic's points taken for quadratic k are those of its
piece between knots t_k and t_{k+1} (the knots the scheme used), at
u = 0, 1/8, ..., 1.

    g_estimate(C, Q, knots, nodes=9) -> (N,)

No nearest-point search, no iteration: polynomial evaluations at fixed nodes.
It is an estimate, not a bound.  It returns infinity for a spline with a
quadratic whose control point projects outside its own chord (such a quadratic
overshoots an end, which the distance to the parabola cannot see).
"""
import numpy as np
from common import power_basis
from schemes import cross


def g_estimate(C, Q, knots, nodes=9):
    C = np.asarray(C, dtype=float); Q = np.asarray(Q, dtype=float)
    N, n = Q.shape[:2]
    on = np.concatenate([C[:, 0:1], (Q[:, :-1] + Q[:, 1:]) / 2, C[:, 3:4]], axis=1)
    b0, b1, b2, b3 = power_basis(C)
    u = np.linspace(0.0, 1.0, nodes)
    t = knots[:, :-1, None] + (knots[:, 1:, None] - knots[:, :-1, None]) * u[None, None, :]     # (N, n, U)
    tt = t[..., None]
    P = b0[:, None, None] + tt * (b1[:, None, None] + tt * (b2[:, None, None] + tt * b3[:, None, None]))
    a = (on[:, :-1] - Q)[:, :, None, :]; b = (on[:, 1:] - Q)[:, :, None, :]
    w = P - Q[:, :, None, :]
    axb = cross(a, b)
    la = np.hypot(a[..., 0], a[..., 1]); lb = np.hypot(b[..., 0], b[..., 1])
    flat = np.abs(axb) <= 1e-9 * (la * lb + 1e-300)
    axb_s = np.where(flat, 1.0, axb)
    x = cross(w, b) / axb_s
    y = cross(a, w) / axb_s
    G = (1 - x - y) ** 2 - 4 * x * y
    Gx = 2 * x - 2 * y - 2; Gy = 2 * y - 2 * x - 2
    v = Gx[..., None] * b - Gy[..., None] * a
    den = np.hypot(v[..., 0], v[..., 1])
    d = np.abs(G) * np.abs(axb_s) / np.where(den > 1e-300, den, 1.0)
    # a flat quadratic is a segment of the line through its ends: distance to that line
    chord = (on[:, 1:] - on[:, :-1])[:, :, None, :]
    cl = np.hypot(chord[..., 0], chord[..., 1])
    dline = np.abs(cross(P - on[:, :-1, None, :], chord)) / np.where(cl > 1e-300, cl, 1.0)
    dpt = np.hypot(*(np.moveaxis(P - on[:, :-1, None, :], -1, 0)))
    # ... but along that line the two can still differ (an overshoot), so for a flat
    # quadratic take the larger of that and the distance at equal parameters
    uu = u[None, None, :, None]
    Sq = (1 - uu) ** 2 * on[:, :-1, None, :] + 2 * uu * (1 - uu) * Q[:, :, None, :] + uu ** 2 * on[:, 1:, None, :]
    dep = np.hypot(*(np.moveaxis(P - Sq, -1, 0)))
    d = np.where(flat, np.maximum(np.where(cl > 1e-12, dline, dpt), dep), d)
    est = d.max(axis=(1, 2))
    # A quadratic whose control point projects outside its chord runs past one of its
    # ends and comes back; the distance to the parabola does not see that.  No estimate then.
    c2 = chord[:, :, 0, :]
    back = ((Q - on[:, :-1]) * c2).sum(-1) < 0
    fwd = ((on[:, 1:] - Q) * c2).sum(-1) < 0
    est = np.where((back | fwd).any(1), np.inf, est)
    return est


if __name__ == '__main__':
    import json, os
    import schemes
    from common import CORPUS
    BIG = '/home/behdad/hb/cu2qu-work/research-weights'
    C = np.array(json.load(open(CORPUS)), dtype=float)
    N = len(C)
    print('node-evaluated estimate against the measured Hausdorff distance, whole corpus (measured > 0.01):')
    print('   scheme  n   cubics   estimate/measured: min     1%      10%     median   90%     99%     max')
    for key, build in (('ramp', lambda n: (schemes.ramp(C, n) if n > 1 else schemes.tangents(C, schemes.uniform_knots(N, 1)), schemes.uniform_knots(N, n))),
                       ('optu', lambda n: (schemes.optu(C, n) if n > 1 else schemes.tangents(C, schemes.uniform_knots(N, 1)), schemes.uniform_knots(N, n))),
                       ('zero', lambda n: (schemes.zero(C, n), schemes.uniform_knots(N, n))),
                       ('Tgam', lambda n: (schemes.tangents(C, schemes.gamma_knots(C, n)), schemes.gamma_knots(C, n)))):
        H = np.load(os.path.join(BIG, f'corpus_{key}.npy'))
        for n in (1, 2, 3, 4):
            Q, kn = build(n)
            est = g_estimate(C, Q, kn)
            h = H[:, n - 1]
            ok = np.isfinite(est) & np.isfinite(h) & (h > 0.01)
            r = est[ok] / h[ok]
            q = np.quantile(r, [0.01, 0.1, 0.5, 0.9, 0.99])
            print(f'   {key:6s} {n:2d}  {ok.sum():6d}   {r.min():6.3f}  ' + '  '.join(f'{v:6.3f}' for v in q) + f'  {r.max():8.3f}')
