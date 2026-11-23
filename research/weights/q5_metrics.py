#!/usr/bin/env python3
"""Q5: side effects of each scheme on the corpus, at the n the scheme picks.

  tangent error   angle between the spline's and the cubic's direction at an end
                  (both ends pooled; a spline whose end off-curve point coincides
                  with the end point is taken to leave along its next point)
  on-curve offset largest distance of an implied on-curve point from the cubic
  extra inflections  sign changes of the turning direction between consecutive
                  quadratics beyond the number of inflections of the cubic
  curvature jump  at the on-curve points of inflection-free cubics: curvature
                  just after / just before = area of triangle k / area of
                  triangle k-1 (a quadratic spline always jumps; the cubic's own
                  curvature also changes from piece to piece, so compare schemes)

d3-family schemes use the n of their equal-parameter formula; the geometric
schemes (T, BT-kept, BT-free) the n of geo2.py (estimate with a 1.1 margin).
"""
import json, os, sys
import numpy as np
import schemes, chd
from common import CORPUS, phi_max, ramp as ramp_w
from schemes import cross, optu_weights
from corpus_stats import formula_n, BIG

C = np.array(json.load(open(CORPUS)), dtype=float)
N = len(C)
d3v = C[:, 3] - 3 * C[:, 2] + 3 * C[:, 1] - C[:, 0]
d3 = np.hypot(d3v[:, 0], d3v[:, 1])
ft_n = np.load(os.path.join(BIG, 'corpus_ft_n.npy'))
e0, e1 = schemes.unit_end_tangents(C)
INFL = schemes.inflections(C)


def point_to_cubic(Cs, P):
    """distance from points P (M, 2) to their cubics: 257 samples of [0, 1], then 129 around the best."""
    b0 = Cs[:, 0]; b1 = 3 * (Cs[:, 1] - Cs[:, 0]); b2 = 3 * (Cs[:, 0] - 2 * Cs[:, 1] + Cs[:, 2]); b3 = Cs[:, 3] - 3 * Cs[:, 2] + 3 * Cs[:, 1] - Cs[:, 0]
    def dist(tt):
        tt = tt[..., None]
        Bp = b0[:, None] + tt * (b1[:, None] + tt * (b2[:, None] + tt * b3[:, None]))
        return np.sqrt(((Bp - P[:, None, :]) ** 2).sum(-1))
    t1 = np.broadcast_to(np.linspace(0, 1, 257), (len(P), 257))
    d1 = dist(t1)
    tb = t1[np.arange(len(P)), d1.argmin(1)]
    best = d1.min(1)
    lo = np.clip(tb - 1 / 256, 0, 1); hi = np.clip(tb + 1 / 256, 0, 1)
    t = tb.copy()
    for _ in range(8):          # Newton on (B(t) - P) . B'(t), kept inside the bracket
        tt = t[:, None]
        B = b0 + tt * (b1 + tt * (b2 + tt * b3)); D = b1 + tt * (2 * b2 + 3 * tt * b3); DD = 2 * b2 + 6 * tt * b3
        f = ((B - P) * D).sum(-1); fp = (D * D).sum(-1) + ((B - P) * DD).sum(-1)
        t = np.clip(t - np.where(fp > 1e-30, f / np.where(fp > 1e-30, fp, 1.0), 0.0), lo, hi)
    tt = t[:, None]
    B = b0 + tt * (b1 + tt * (b2 + tt * b3))
    return np.minimum(best, np.sqrt(((B - P) ** 2).sum(-1)))


def metrics(groups):
    """groups: list of (sel, Q) with Q (len(sel), n, 2), n >= 1."""
    ang = []; off = []; extra = 0; with_extra = 0; jumps = []; total = 0
    for sel, Q in groups:
        n = Q.shape[1]
        total += len(sel)
        Cs = C[sel]
        on = np.concatenate([Cs[:, 0:1], (Q[:, :-1] + Q[:, 1:]) / 2, Cs[:, 3:4]], axis=1)
        # tangent error at both ends
        for P, Qe, nxt, t in ((Cs[:, 0], Q[:, 0], on[:, 1], e0[sel]), (Cs[:, 3], Q[:, -1], on[:, -2], -e1[sel])):
            v = Qe - P
            small = np.hypot(v[:, 0], v[:, 1]) < 1e-9
            v = np.where(small[:, None], nxt - P, v)
            ok = np.hypot(v[:, 0], v[:, 1]) > 1e-9
            a = np.degrees(np.arctan2(np.abs(cross(t, v)), (t * v).sum(-1)))
            ang.append(a[ok])
        if n >= 2:
            best = np.zeros(len(sel))
            for k in range(1, n):
                best = np.maximum(best, point_to_cubic(Cs, on[:, k]))
            off.append(best)
            tr = cross(Q - on[:, :-1], on[:, 1:] - Q)
            scale = np.abs(tr).max(1, keepdims=True) + 1e-300
            sg = np.sign(np.where(np.abs(tr) > 1e-9 * scale, tr, 0.0))
            changes = ((sg[:, 1:] * sg[:, :-1]) < 0).sum(1)
            ex = np.maximum(changes - INFL[sel], 0)
            extra += int(ex.sum()); with_extra += int((ex > 0).sum())
            okc = (INFL[sel] == 0) & (np.abs(tr) > 1e-9 * scale).all(1)
            if okc.any():
                ratio = np.abs(tr[okc][:, 1:] / tr[okc][:, :-1])
                jumps.append(np.abs(np.log(ratio)).max(1))
        else:
            off.append(np.zeros(len(sel)))
    ang = np.concatenate(ang); off = np.concatenate(off)
    jumps = np.concatenate(jumps) if jumps else np.zeros(1)
    q = np.quantile
    return (f'{q(ang, .5):6.3f} {q(ang, .9):6.3f} {q(ang, .99):6.3f} {ang.max():7.2f} | '
            f'{q(off, .5):5.3f} {q(off, .9):5.3f} {off.max():6.3f} | {extra:4d} in {with_extra:4d} | '
            f'{np.exp(q(jumps, .5)):5.2f} {np.exp(q(jumps, .9)):5.2f} {np.exp(q(jumps, .99)):6.2f}')


def d3_groups(weights_of, n_of):
    out = []
    for n in np.unique(n_of):
        sel = np.where(n_of == n)[0]
        if n == 1:      # tangents kept: cu2qu's single quadratic; otherwise the family's
            w = weights_of(1)
            Q = schemes.tangents(C[sel], schemes.uniform_knots(len(sel), 1)) if w is None else schemes.family_offcurves(C[sel], w)
        else:
            Q = schemes.family_offcurves(C[sel], weights_of(int(n)))
        out.append((sel, Q))
    return out


if __name__ == '__main__':
    print(f'{N} cubics; with 0 / 1 / 2 inflections: {(INFL == 0).sum()} / {(INFL == 1).sum()} / {(INFL == 2).sum()}')
    Phi_ramp = lambda n: phi_max(ramp_w(n))
    Phi_opt = lambda n: 4 / 27 if n % 2 == 0 else phi_max(optu_weights(n))
    Phi_zero = lambda n: 3 ** 0.5 / 36
    Phi_r50 = lambda n: phi_max(optu_weights(n, 0.5)) if n > 1 else 3 ** 0.5 / 36
    import geo2
    hdr = ('scheme                      tol | end tangent error, degrees     | on-curve offset      | extra        | curvature jump factor\n'
           '                                |  median    90%    99%     max  | median  90%    max   | inflections  | median  90%    99%')
    print(hdr)
    for j, tol in enumerate((1.0, 0.5, 0.1)):
        one = ft_n[:, j] == 1
        n_ramp = np.where(one, 1, formula_n(d3, tol, Phi_ramp, start=2))
        n_opt = np.where(one, 1, formula_n(d3, tol, Phi_opt, start=2))
        rows = [
            ('ramp (cu2qu)', d3_groups(lambda n: None if n == 1 else ramp_w(n), n_ramp)),
            ('optu, its own n', d3_groups(lambda n: None if n == 1 else optu_weights(n), n_opt)),
            ('optu at the ramp\'s n', d3_groups(lambda n: None if n == 1 else optu_weights(n), n_ramp)),
            ('relax50', d3_groups(lambda n: optu_weights(n, 0.5), formula_n(d3, tol, Phi_r50))),
            ('zero', d3_groups(lambda n: np.zeros(n), formula_n(d3, tol, Phi_zero))),
        ]
        for kind in ('T', 'BT-kept', 'BT-free'):
            n_out, by_est, H, groups = geo2.run(kind, tol, j, want_groups=True)
            rows.append((kind + ' (geo2)', groups))
        for name, groups in rows:
            print(f'{name:26s} {tol:4.1f} | {metrics(groups)}', flush=True)
        print()
