#!/usr/bin/env python3
"""Generate random and adversarial cubics (N x 8 doubles each) into W/sets/.

  random      four control points uniform in a square of side 10^U(0.5,3.5)
  fontlike    chord 20..600 units, handles within +-70 degrees of the chord, lengths 0.05..0.9 of it, integer coordinates
  arcs        circular-arc-like cubics of 10..170 degrees with perturbed handle lengths
  cusp        cubics with an exact or near cusp (B'(t0) = 0 for some t0 in (0,1)), perturbed
  loop        self-intersecting cubics (crossed handles)
  inflect     S-shaped cubics, including ones whose end tangents are nearly parallel
  coincident  c1 = c0, c2 = c3, c1 = c2, three equal, all equal, c0 = c3, and near-coincident versions
  straight    collinear and nearly collinear control points with uneven handles, some overshooting the ends
  parallel    end tangents nearly parallel (tangent intersection far away)
  tiny_huge   the 'fontlike' shapes scaled by 10^U(-3, 4)
"""
import os, sys
import numpy as np
from common import W

N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261122
out = W + (sys.argv[3] if len(sys.argv) > 3 else '/sets')     # usage: gen_sets.py [N [seed [/subdir]]]
rng = np.random.default_rng(seed)
os.makedirs(out, exist_ok=True)

def rot(v, ang):
    c, s = np.cos(ang), np.sin(ang)
    return np.stack([c * v[..., 0] - s * v[..., 1], s * v[..., 0] + c * v[..., 1]], -1)

def place(c):
    """random rotation + translation"""
    ang = rng.uniform(0, 2 * np.pi, len(c))[:, None]          # (n, 1) against (n, 4)
    c = rot(c, ang)
    assert c.shape == (len(ang), 4, 2)
    return c + rng.uniform(-500, 500, (len(c), 1, 2))

def save(name, c):
    c = np.asarray(c, float).reshape(-1, 4, 2)
    c.reshape(-1, 8).tofile(f'{out}/{name}.bin')
    print(name, len(c))

def gen_random(n):
    S = 10 ** rng.uniform(0.5, 3.5, (n, 1, 1))
    return rng.uniform(-1, 1, (n, 4, 2)) * S

def gen_fontlike(n, integer=True):
    L = rng.uniform(20, 600, n)
    a0 = rng.uniform(-70, 70, n) * np.pi / 180
    a1 = rng.uniform(-70, 70, n) * np.pi / 180
    # mostly C-shaped: same-side handles
    same = rng.random(n) < 0.8
    a1 = np.where(same, -np.sign(a0) * np.abs(a1), a1)
    l0 = rng.uniform(0.05, 0.9, n) * L
    l1 = rng.uniform(0.05, 0.9, n) * L
    c = np.zeros((n, 4, 2))
    c[:, 1] = np.stack([l0 * np.cos(a0), l0 * np.sin(a0)], -1)
    c[:, 3, 0] = L
    c[:, 2] = c[:, 3] - np.stack([l1 * np.cos(a1), -l1 * np.sin(a1)], -1) * np.array([1, -1])
    c = place(c)
    return np.round(c) if integer else c

def gen_arcs(n):
    R = 10 ** rng.uniform(1, 3, n)
    th = rng.uniform(10, 170, n) * np.pi / 180
    k = 4 / 3 * np.tan(th / 4)
    k0 = k * rng.uniform(0.6, 1.5, n); k1 = k * rng.uniform(0.6, 1.5, n)
    c = np.zeros((n, 4, 2))
    c[:, 0] = np.stack([R, 0 * R], -1)
    c[:, 1] = np.stack([R, k0 * R], -1)
    e = np.stack([R * np.cos(th), R * np.sin(th)], -1)
    te = np.stack([-np.sin(th), np.cos(th)], -1)
    c[:, 3] = e
    c[:, 2] = e - te * (k1 * R)[:, None]
    return place(c)

def gen_cusp(n):
    # B'(t0) = 0: choose b1 + 2 b2 t0 + 3 b3 t0^2 = 0 with random b2, b3, t0
    t0 = rng.uniform(0.05, 0.95, n)[:, None]
    b2 = rng.uniform(-300, 300, (n, 2)); b3 = rng.uniform(-300, 300, (n, 2))
    b1 = -(2 * b2 * t0 + 3 * b3 * t0 * t0)
    # perturb: exact cusp, near cusp (small loop or small bump)
    pert = np.where(rng.random((n, 1)) < 0.3, 0.0, 10 ** rng.uniform(-6, 0.5, (n, 1)))
    b1 = b1 + pert * rng.normal(size=(n, 2))
    c0 = rng.uniform(-200, 200, (n, 2))
    c1 = c0 + b1 / 3; c2 = c1 + (b2 + b1) / 3; c3 = c0 + b1 + b2 + b3
    return np.stack([c0, c1, c2, c3], 1)

def gen_loop(n):
    L = rng.uniform(20, 400, n)
    c = np.zeros((n, 4, 2))
    c[:, 3, 0] = L
    c[:, 1] = np.stack([L * rng.uniform(0.8, 2.5, n), L * rng.uniform(0.2, 2.0, n)], -1)
    c[:, 2] = np.stack([L * rng.uniform(-1.5, 0.2, n), L * rng.uniform(0.2, 2.0, n)], -1)
    return place(c)

def gen_inflect(n):
    L = rng.uniform(20, 600, n)
    a0 = rng.uniform(1, 70, n) * np.pi / 180
    par = rng.random(n) < 0.4                     # nearly parallel end tangents
    a1 = np.where(par, a0 + 10 ** rng.uniform(-6, -1, n) * rng.choice([-1, 1], n), rng.uniform(1, 70, n) * np.pi / 180)
    l0 = rng.uniform(0.05, 0.9, n) * L; l1 = rng.uniform(0.05, 0.9, n) * L
    c = np.zeros((n, 4, 2))
    c[:, 1] = np.stack([l0 * np.cos(a0), l0 * np.sin(a0)], -1)
    c[:, 3, 0] = L
    c[:, 2] = c[:, 3] - np.stack([l1 * np.cos(a1), l1 * np.sin(a1)], -1)
    return place(c)

def gen_coincident(n):
    base = gen_fontlike(n, integer=False)
    kind = rng.integers(0, 9, n)
    c = base.copy()
    eps = np.where(rng.random((n, 1)) < 0.5, 0.0, 10 ** rng.uniform(-9, -1, (n, 1))) * rng.normal(size=(n, 2))
    for i in range(n):
        k = kind[i]
        if k == 0: c[i, 1] = c[i, 0] + eps[i]
        elif k == 1: c[i, 2] = c[i, 3] + eps[i]
        elif k == 2: c[i, 1] = c[i, 0] + eps[i]; c[i, 2] = c[i, 3] - eps[i]
        elif k == 3: c[i, 2] = c[i, 1] + eps[i]
        elif k == 4: c[i, 1] = c[i, 0]; c[i, 2] = c[i, 0] + eps[i]
        elif k == 5: c[i, 1] = c[i, 3] + eps[i]; c[i, 2] = c[i, 3]
        elif k == 6: c[i, 1] = c[i, 0]; c[i, 2] = c[i, 0]; c[i, 3] = c[i, 0] + eps[i]
        elif k == 7: c[i, 3] = c[i, 0] + eps[i]
        elif k == 8: c[i, 1] = c[i, 2] = (c[i, 0] + c[i, 3]) / 2 + eps[i]
    return c

def gen_straight(n):
    L = rng.uniform(5, 800, n)
    t1 = rng.choice([0.001, 0.01, 0.05, 0.2, 0.5, 0.9, 1.2, -0.2], n) * rng.uniform(0.5, 1.5, n)
    t2 = rng.choice([0.999, 0.99, 0.95, 0.8, 0.5, 0.1, -0.2, 1.3], n) * rng.uniform(0.5, 1.5, n)
    noise = np.where(rng.random((n, 1)) < 0.25, 0.0, 10 ** rng.uniform(-6, 0.7, (n, 1)))
    c = np.zeros((n, 4, 2))
    c[:, 1, 0] = t1 * L; c[:, 2, 0] = t2 * L; c[:, 3, 0] = L
    c[:, 1, 1] = noise[:, 0] * rng.normal(size=n); c[:, 2, 1] = noise[:, 0] * rng.normal(size=n)
    c = place(c)
    rnd = rng.random(n) < 0.3
    c[rnd] = np.round(c[rnd])
    return c

def gen_parallel(n):
    L = rng.uniform(20, 600, n)
    a0 = rng.uniform(-80, 80, n) * np.pi / 180
    da = 10 ** rng.uniform(-8, -0.5, n) * rng.choice([-1, 1], n)
    l0 = rng.uniform(0.05, 0.9, n) * L; l1 = rng.uniform(0.05, 0.9, n) * L
    c = np.zeros((n, 4, 2))
    c[:, 1] = np.stack([l0 * np.cos(a0), l0 * np.sin(a0)], -1)
    c[:, 3, 0] = L
    # end tangent direction = start tangent direction rotated by da (U-shape: reversed)
    rev = rng.random(n) < 0.5
    a1 = np.where(rev, a0 + np.pi + da, a0 + da)
    c[:, 2] = c[:, 3] - np.stack([l1 * np.cos(a1), l1 * np.sin(a1)], -1)
    return place(c)

def gen_tiny_huge(n):
    c = gen_fontlike(n, integer=False)
    s = 10 ** rng.uniform(-3, 4, (n, 1, 1))
    return c * s

save('random', gen_random(N))
save('fontlike', gen_fontlike(N))
save('arcs', gen_arcs(N))
save('cusp', gen_cusp(N))
save('loop', gen_loop(N))
save('inflect', gen_inflect(N))
save('coincident', gen_coincident(N))
save('straight', gen_straight(N))
save('parallel', gen_parallel(N))
save('tiny_huge', gen_tiny_huge(N))
