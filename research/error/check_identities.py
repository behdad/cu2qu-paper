#!/usr/bin/env python3
"""Numerical checks of the identities the proofs rest on, on random cubics.

 1. Quadratic k minus cubic piece k at equal local parameter u is d phi_k(u), d = d3 / n^3
    (ramp, s = 0, alternating), and the spline is C1 with legs  a_k = B'(t_k)/(2n) - d (2 g_k + h_k)/2.
 2. Lemma N (normal line): for p = Q(u) + E, with w = V/|V|, eps = (E.w)/(|V| (1 + sqrt(1 + y))),
    y = ((z - a).w)(E.w)/|V|^2:   Q(u + eps) - p = -(E.n - (C/|V|) eps^2) n.
 3. Lemma G (gradient line): G = xi^2 - eta vanishes on the parabola, grad G = -rot(V(xi))/C, and the
    point p - sign(G) s g/|g| with s = c |G| |C| / |V|, c = 2 / (1 + sqrt(1 - G (V.W)^2/|V|^4)), is on it.
"""
import numpy as np
from geo import *

rng = np.random.default_rng(11)
N = 4000
c = rng.uniform(-500, 500, (N, 4, 2))
scale = 500.0

# 1
worst = 0; worst_leg = 0
for mode in ('ramp', 'zero', 'alt'):
    for n in (2, 3, 4, 5, 8, 13):
        if mode == 'alt' and n % 2:
            continue
        s = weights(mode, n)
        off, on, d = spline(c, s)
        co, g = phi_coeffs(s)
        h = (1 - s) / 2
        c0, c1, c2, c3 = c[:, 0], c[:, 1], c[:, 2], c[:, 3]
        b1 = 3 * (c1 - c0); b2 = 3 * (c0 - 2 * c1 + c2); b3 = c3 - 3 * c2 + 3 * c1 - c0
        for k in range(n):
            tk = k / n
            Bp = b1 + tk * (2 * b2 + 3 * tk * b3)
            a_pred = Bp / (2 * n) - d * (2 * g[k] + h[k]) / 2
            worst_leg = max(worst_leg, np.abs(off[:, k] - on[:, k] - a_pred).max() / scale)
            for u in np.linspace(0, 1, 9):
                t = (k + u) / n
                B = c0 + t * (b1 + t * (b2 + t * b3))
                Q = (1 - u) ** 2 * on[:, k] + 2 * u * (1 - u) * off[:, k] + u * u * on[:, k + 1]
                phi = co[k, 0] + u * (co[k, 1] + u * (co[k, 2] + u * co[k, 3]))
                worst = max(worst, np.abs(Q - B - d * phi).max() / scale)
print(f'1. quadratic - cubic = d phi_k(u): largest residual / size {worst:.2e};  leg formula: {worst_leg:.2e}')

# 2 and 3: random parabolas and random points near them
P0 = rng.uniform(-100, 100, (N, 2)); a = rng.uniform(-100, 100, (N, 2)); z = rng.uniform(-100, 100, (N, 2))
u = rng.uniform(0, 1, N)
E = rng.normal(size=(N, 2)) * 10 ** rng.uniform(-3, 0.5, (N, 1))
V = (1 - u)[:, None] * a + u[:, None] * z
Vn = np.sqrt(dot(V, V)); w = V / Vn[:, None]; nrm = np.stack([-w[:, 1], w[:, 0]], -1)
C = cross(a, z)
Q = lambda x: P0 + 2 * x[:, None] * a + (x * x)[:, None] * (z - a)
p = Q(u) + E
Ew = dot(E, w); En = dot(E, nrm)
y = dot(z - a, w) * Ew / Vn ** 2
ok = y > -1
eps = Ew / (Vn * (1 + np.sqrt(np.where(ok, 1 + y, 1))))
res = Q(u + eps) - p + (En - C / Vn * eps ** 2)[:, None] * nrm
print(f'2. normal-line lemma: {ok.sum()} cases with y > -1, largest residual {np.abs(res[ok]).max():.2e}')

X = cross(p - P0, z) / C; Y = cross(a, p - P0) / C          # p = P0 + X a + Y z
xi = (X + Y) / 2; G = xi * xi - Y
Vx = (1 - xi)[:, None] * a + xi[:, None] * z
Vxn = np.sqrt(dot(Vx, Vx))
grad = -np.stack([-Vx[:, 1], Vx[:, 0]], -1) / C[:, None]
# finite-difference check of the gradient
hh = 1e-5
def Gf(q):
    Xq = cross(q - P0, z) / C; Yq = cross(a, q - P0) / C
    return ((Xq + Yq) / 2) ** 2 - Yq
gnum = np.stack([(Gf(p + [hh, 0]) - Gf(p - [hh, 0])) / (2 * hh), (Gf(p + [0, hh]) - Gf(p - [0, hh])) / (2 * hh)], -1)
print(f'3. gradient of G: largest relative error vs finite differences {(np.abs(gnum - grad).max(1) / np.abs(grad).max(1)).max():.2e}')
Wv = a - z
chi = G * dot(Vx, Wv) ** 2 / Vxn ** 4
ok = chi < 1
cf = 2 / (1 + np.sqrt(np.where(ok, 1 - chi, 1)))
sdist = cf * np.abs(G) * np.abs(C) / Vxn
ghat = grad / np.sqrt(dot(grad, grad))[:, None]
hit = p - (np.sign(G) * sdist)[:, None] * ghat
print(f'   gradient-line lemma: {ok.sum()} cases with chi < 1, largest |G| at the hit point relative to |G(p)|: {np.abs(Gf(hit)[ok] / np.where(G[ok] == 0, 1, G[ok])).max():.2e}')
xi_hit = (cross(hit - P0, z) / C + cross(a, hit - P0) / C) / 2
dxi = cf * np.abs(G) * np.abs(dot(Vx, Wv)) / (2 * Vxn ** 2)
print(f'   parameter of the hit point: |xi_hit - xi| against c |G| |V.W| / (2 |V|^2): largest difference {np.abs(np.abs(xi_hit - xi) - dxi)[ok].max():.2e}')
