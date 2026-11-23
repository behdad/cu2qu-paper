#!/usr/bin/env python3
"""Balanced tangent-intersection splines ("BT"): closed-form constructions that
follow the fourth-order model.

    bt(C, knots, tangent)   off-curve points (N, n, 2), validity (N,), estimate (N,)

Start from the T spline at the given knots.  For piece k let e_k be the signed
first-order error of its T quadratic (tquad.t_profile, the value of largest
size; positive = the quadratic is on the outside, the side of Q), E0 = 16 e_k,
theta_k the angle the piece turns by.  Move off-curve point k by

    beta_k E0  toward the curve (along the normal of the chord), and
    - alpha_k E0 / theta_k  forward along the chord,

with the universal beta_k, alpha_k of q2_balance_weights.py (they depend on n,
k and on whether the end tangents are kept).  With tangents kept the two end
points slide along the end tangents instead (by - alpha E0 / theta).

Estimate of the result's error: c(n) * 16 * max_k |e_k|, c(n) the model
constant.  A cubic with a piece that has no valid tangent intersection is
flagged invalid (the caller falls back to a d3-family scheme).
"""
import numpy as np
import schemes, tquad
from schemes import cross
from q2_balance_weights import solve as model_solve

_TAB = {}


def table(n, tangent, variant='full'):
    key = (n, tangent, variant)
    if key not in _TAB:
        _TAB[key] = model_solve(n, tangent, variant)
    return _TAB[key]


def bt(C, knots, tangent, variant='full', U=17):
    C = np.asarray(C, dtype=float)
    N = len(C); n = knots.shape[1] - 1
    c, beta, alpha = table(n, tangent, variant)
    pc = tquad.t_pieces(C, knots)
    D = tquad.t_profile(pc, np.ones((N, n)), np.ones((N, n)), U)       # each piece as its own G1 quadratic
    j = np.abs(D).argmax(-1)
    e = np.take_along_axis(D, j[..., None], -1)[..., 0]               # signed, (N, n)
    valid = pc['valid'].all(1)
    A = pc['A']; B = pc['B']; Q = pc['Q']
    chord = B - A
    cl = np.hypot(chord[..., 0], chord[..., 1]); cl = np.where(cl > 1e-12, cl, 1.0)
    That = chord / cl[..., None]
    # unit vector from Q toward the chord line
    w = A - Q
    perp = w - (w * That).sum(-1)[..., None] * That
    h = np.hypot(perp[..., 0], perp[..., 1])
    nhat = perp / np.where(h > 1e-12, h, 1.0)[..., None]
    a = pc['a']; b = pc['b']
    la = np.hypot(a[..., 0], a[..., 1]); lb = np.hypot(b[..., 0], b[..., 1])
    sinth = np.abs(cross(a, b)) / np.where(la * lb > 1e-12, la * lb, 1.0)       # sin of the turning angle
    cosm = -(a * b).sum(-1) / np.where(la * lb > 1e-12, la * lb, 1.0)            # cos of the turning angle
    theta = np.arctan2(sinth, cosm)
    E0 = 16 * e
    Qn = Q + (beta[None, :] * E0)[..., None] * nhat
    tang = np.where(theta > 1e-9, -alpha[None, :] * E0 / np.where(theta > 1e-9, theta, 1.0), 0.0)
    Qn = Qn + tang[..., None] * That
    if tangent:
        e0, e1 = schemes.unit_end_tangents(C)
        Qn[:, 0] = Q[:, 0] + tang[:, 0, None] * e0
        if n > 1:
            Qn[:, -1] = Q[:, -1] + tang[:, -1, None] * e1
        else:
            Qn[:, 0] = Q[:, 0]
    est = c * 16 * np.abs(e).max(1)
    return Qn, valid, est


def single(C, rule='N'):
    """One quadratic, tangents free.
    'N':  Q* = Q_T - 2 (sqrt 2 - 1) P_N (Q_T - q_M),  P_N the projection on the cubic's normal at t = 1/2
    'D':  bt() with n = 1"""
    C = np.asarray(C, dtype=float)
    if rule == 'D':
        return bt(C, schemes.uniform_knots(len(C), 1), False)
    QT, valid = schemes.tangents(C, schemes.uniform_knots(len(C), 1), return_valid=True)
    QT = QT[:, 0]; valid = valid[:, 0]
    qM = (3 * (C[:, 1] + C[:, 2]) - C[:, 0] - C[:, 3]) / 4
    T = C[:, 3] + C[:, 2] - C[:, 1] - C[:, 0]
    L = np.hypot(T[:, 0], T[:, 1])
    Nrm = np.stack([-T[:, 1], T[:, 0]], axis=1) / np.where(L > 1e-12, L, 1.0)[:, None]
    comp = ((QT - qM) * Nrm).sum(-1)
    Qs = QT - 2 * (2 ** 0.5 - 1) * comp[:, None] * Nrm
    est = (3 - 2 * 2 ** 0.5) * np.abs(comp) / 2
    return Qs[:, None, :], valid, est


if __name__ == '__main__':
    import math
    import chd
    from q2_model_check import arc_cubic
    from q2_fourth_order import solve as solve_c1
    print('closed-form balanced T splines on circular arcs: measured error / E0 against the model constant')
    for deg in (20, 60):
        Theta = math.radians(deg); R = 1000.0
        c = arc_cubic(R, Theta); C = c[None]
        print(f'  arc of {deg} degrees')
        for tangent in (True, False):
            row = []
            for n in range(1, 7):
                E0 = R * (Theta / n) ** 4 / 8
                Qn, valid, est = bt(C, schemes.affine_knots(C, n), tangent)
                h = chd.hd1(c, Qn[0], m=64)[0]
                ref, _ = solve_c1(n, tangent)
                row.append(f'{h / E0:.4f}/{ref:.4f} (est {est[0] / E0:.4f})')
            print('    ' + ('kept: ' if tangent else 'free: ') + '  '.join(row))
