#!/usr/bin/env python3
"""n = 1: the cubic against its tangent-intersection quadratic.

Checks, on the corpus and on random cubics:

 1. The identity.  With q the intersection of the end tangents, a = q - c0,
    z = c3 - q, c1 = c0 + l1 a, c2 = c3 - l2 z, p = (3 l1 - 2)/6,
    r = (2 - 3 l2)/6, e(t) = (1-t) p + t r:

        B(t) = Q(tau(t)) + rho(t) (a - z)
        tau(t) = t + 3 t (1-t) e(t)
        rho(t) = 3 t^2 (1-t)^2 (2 (p - r) + 3 e(t)^2)

    rho is the implicit polynomial G of the parabola evaluated on the cubic.

 2. The level-set ("Sampson") distance |G| / |grad G| = |rho| |a x z| / |V(tau)|,
    V(tau) = (1-tau) a + tau z, with the factor c for points outside the
    parabola, maximised over t by dense sampling, against the true
    Hausdorff distance measured by gt.c.
"""
import numpy as np
from common import *

def cross(u, v):
    return u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]

def dot(u, v):
    return (u * v).sum(-1)

def n1_setup(c):
    """fontTools' tangent intersection, and the quantities of the identity."""
    c0, c1, c2, c3 = c[:, 0], c[:, 1], c[:, 2], c[:, 3]
    ab = c1 - c0; cd = c3 - c2
    pp = np.stack([-ab[:, 1], ab[:, 0]], -1)
    den = dot(pp, cd)
    num = dot(pp, c0 - c2)
    with np.errstate(divide='ignore', invalid='ignore'):
        h = np.where(np.abs(den) < 1e-15, np.nan, num / den)
    q = c2 + cd * h[:, None]
    a = q - c0; z = c3 - q
    with np.errstate(divide='ignore', invalid='ignore'):
        l1 = dot(c1 - c0, a) / dot(a, a)
        l2 = dot(c3 - c2, z) / dot(z, z)
    p = (3 * l1 - 2) / 6
    r = (2 - 3 * l2) / 6
    return q, a, z, l1, l2, p, r

def sampson_dense(c, T=2001):
    q, a, z, l1, l2, p, r = n1_setup(c)
    C = cross(a, z)
    W = a - z
    t = np.linspace(0, 1, T)[None, :]
    e = (1 - t) * p[:, None] + t * r[:, None]
    tau = t + 3 * t * (1 - t) * e
    rho = 3 * t * t * (1 - t) ** 2 * (2 * (p - r)[:, None] + 3 * e * e)
    V = (1 - tau)[..., None] * a[:, None, :] + tau[..., None] * z[:, None, :]
    Vn = np.sqrt(dot(V, V))
    first = np.abs(rho) * np.abs(C)[:, None] / Vn
    x = rho * dot(V, W[:, None, :]) ** 2 / Vn ** 4
    with np.errstate(invalid='ignore'):
        cfac = np.where(rho > 0, 2 / (1 + np.sqrt(np.clip(1 - x, 0, None))), 1.0)
    bad_disc = (x > 1).any(1)
    return (first * cfac).max(1), first.max(1), bad_disc, (l1, l2, p, r)

if __name__ == '__main__':
    c = load_cubics()
    g = load_gt(W + '/gt_ramp.bin', 16)
    H1 = np.fmax(g[:, 0, 1], g[:, 0, 2])
    ep1 = g[:, 0, 0]

    # 1. identity
    q, a, z, l1, l2, p, r = n1_setup(c)
    ok = np.isfinite(q).all(1) & np.isfinite(l1) & np.isfinite(l2)
    t = np.linspace(0, 1, 101)[None, :, None]
    cc = c[ok]
    B = (1-t)**3*cc[:, None, 0] + 3*(1-t)**2*t*cc[:, None, 1] + 3*(1-t)*t*t*cc[:, None, 2] + t**3*cc[:, None, 3]
    e = (1 - t) * p[ok][:, None, None] + t * r[ok][:, None, None]
    tau = t + 3 * t * (1 - t) * e
    rho = 3 * t * t * (1 - t) ** 2 * (2 * (p - r)[ok][:, None, None] + 3 * e * e)
    Q = (1-tau)**2*cc[:, None, 0] + 2*(1-tau)*tau*q[ok][:, None, :] + tau*tau*cc[:, None, 3]
    res = B - Q - rho * (a - z)[ok][:, None, :]
    scale = np.abs(cc).max((1, 2)) + np.abs(q[ok]).max(1)
    print('identity B(t) = Q(tau) + rho (a - z): largest residual / scale over the corpus:',
          (np.abs(res).max((1, 2)) / scale).max())

    # 2. level-set distance vs truth
    sd, first, bad, _ = sampson_dense(c)
    valid = ok & (l1 > 0) & (l1 < 2) & (l2 > 0) & (l2 < 2) & ~bad & np.isfinite(H1)
    print('cubics with a tangent intersection:', ok.sum(), ' with handles in (0, 2) and discriminant ok:', valid.sum())
    m = valid & (H1 > 1e-6)
    for name, v in (('level-set with c', sd), ('level-set without c', first)):
        ratio = v[m] / H1[m]
        print(f'{name:22s} / H(1): min {ratio.min():.6f}  percentiles 1/5/25/50/75/95/99: ' +
              ' '.join(f'{x:.4f}' for x in np.percentile(ratio, [1, 5, 25, 50, 75, 95, 99])) + f'  max {ratio.max():.3f}')
        for tol in (1.0, 0.5, 0.1):
            sel = m & (H1 <= 2 * tol) & (H1 >= tol / 4)
            rr = v[sel] / H1[sel]
            print(f'     H(1) in [{tol/4}, {2*tol}]: n={sel.sum()}  min {rr.min():.6f}  median {np.median(rr):.4f}  99% {np.percentile(rr, 99):.4f}  max {rr.max():.4f}')
    for tol in (1.0, 0.5, 0.1):
        print(f'tol {tol}: single quadratic fits: equal-parameter {np.mean(ep1 <= tol):.4f}  true {np.mean(H1 <= tol):.4f}  '
              f'level-set(dense, valid only) {np.mean(valid & (sd <= tol)):.4f}  under-calls {(valid & (sd <= tol) & (H1 > tol)).sum()}')
