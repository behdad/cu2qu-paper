#!/usr/bin/env python3
"""Check the single-quadratic geometric bound:

  1. the identity B(t) = Q(tau(t)) + rho(t) (a - z), exactly, in rational
     arithmetic on a grid large enough for its degree (so: as polynomials);
  2. bound 1, (3/16) max|g2| |a - z|, against a measured two-sided Hausdorff
     distance (point-to-polyline, both ways) on corpus cubics and on
     constructed ones covering the whole range of handle ratios [0, 2],
     zero-length handles, far tangent intersections and extreme scales.

Usage: n1_validate.py cubics.json
"""
import json, math, random, sys
from fractions import Fraction as F
import numpy as np

# ---------------------------------------------------------------- 1. identity
def identity_exact():
    rnd = random.Random(11)
    worst = 0
    for trial in range(3):
        pts = [(F(rnd.randint(-50, 50)), F(rnd.randint(-50, 50))) for _ in range(3)]
        c0, q, c3 = pts
        a = (q[0] - c0[0], q[1] - c0[1])
        z = (c3[0] - q[0], c3[1] - q[1])
        for l1 in (F(0), F(1, 3), F(2, 3), F(5, 4), F(2)):
            for l2 in (F(0), F(1, 2), F(2, 3), F(7, 5), F(2)):
                c1 = (c0[0] + l1 * a[0], c0[1] + l1 * a[1])
                c2 = (c3[0] - l2 * z[0], c3[1] - l2 * z[1])
                p, r = l1 / 2 - F(1, 3), F(1, 3) - l2 / 2
                for k in range(9):  # degree 6 in t: 7 points suffice
                    t = F(k, 8)
                    e = (1 - t) * p + t * r
                    tau = t + 3 * t * (1 - t) * e
                    rho = 3 * t * t * (1 - t) ** 2 * (2 * (p - r) + 3 * e * e)
                    for d in (0, 1):
                        B = (1 - t) ** 3 * c0[d] + 3 * t * (1 - t) ** 2 * c1[d] + 3 * t * t * (1 - t) * c2[d] + t**3 * c3[d]
                        Q = (1 - tau) ** 2 * c0[d] + 2 * tau * (1 - tau) * q[d] + tau * tau * c3[d]
                        diff = B - Q - rho * (a[d] - z[d])
                        assert diff == 0, (trial, l1, l2, t, diff)
                    # tau is the Bezier function with control values 0, l1/2, 1 - l2/2, 1
                    tau2 = 3 * t * (1 - t) ** 2 * (l1 / 2) + 3 * t * t * (1 - t) * (1 - l2 / 2) + t**3
                    assert tau == tau2
    print("identity: exact in rational arithmetic on 3 x 5 x 5 x 9 points (degree 6 in t, 2 in each handle ratio)")


identity_exact()

# ------------------------------------------------------------ 2. the bound
def bound1(c):
    c0, c1, c2, c3 = c
    ab = c1 - c0
    cd = c3 - c2
    if not ab.any():
        ab = c2 - c0          # zero-length handle: the curve leaves toward c2
    if not cd.any():
        cd = c3 - c1
    den = ab[0] * cd[1] - ab[1] * cd[0]
    if den == 0:
        return None, None
    h = (ab[0] * (c0 - c2)[1] - ab[1] * (c0 - c2)[0]) / den
    q = c2 + cd * h
    a = q - c0
    z = c3 - q
    aa, zz = a @ a, z @ z
    if not (aa > 0 and zz > 0):
        return None, None
    l1 = ((c1 - c0) @ a) / aa
    l2 = ((c3 - c2) @ z) / zz
    if not (0 <= l1 <= 2 and 0 <= l2 <= 2):
        return None, None
    p, r = l1 / 2 - 1 / 3, 1 / 3 - l2 / 2
    e2max = max(p * p, r * r)
    e2min = 0.0 if p * r <= 0 else min(p * p, r * r)
    M = max(abs(2 * (p - r) + 3 * e2max), abs(2 * (p - r) + 3 * e2min))
    return 3 / 16 * M * math.hypot(*(a - z)), q


def polyline_dist(P, S):
    """Largest distance from the points P to the polyline S."""
    A, B = S[:-1], S[1:]
    AB = B - A
    L2 = (AB * AB).sum(axis=1)
    L2 = np.where(L2 > 0, L2, 1.0)
    AP = P[:, None, :] - A[None, :, :]
    u = np.clip((AP * AB[None, :, :]).sum(axis=2) / L2[None, :], 0, 1)
    D = AP - u[:, :, None] * AB[None, :, :]
    return np.sqrt((D * D).sum(axis=2)).min(axis=1).max()


K = 600
T = np.linspace(0, 1, K + 1)


def measured(c, q):
    c0, c1, c2, c3 = c
    t = T[:, None]
    B = (1 - t) ** 3 * c0 + 3 * t * (1 - t) ** 2 * c1 + 3 * t * t * (1 - t) * c2 + t**3 * c3
    Q = (1 - t) ** 2 * c0 + 2 * t * (1 - t) * q + t * t * c3
    return max(polyline_dist(B, Q), polyline_dist(Q, B))


def run(name, cubics):
    n = under = 0
    ratios = []
    worst = None
    for c in cubics:
        b, q = bound1(c)
        if b is None:
            continue
        m = measured(c, q)
        size = max(np.abs(c - c[0]).max(), 1e-30)
        n += 1
        # the instrument's own error: the sag of a polyline of K pieces
        if b < m * (1 - 1e-6) - 3e-6 * size:
            under += 1
            if worst is None or m - b > worst[0]:
                worst = (m - b, b, m, c.tolist())
        if m > 1e-5 * size:
            ratios.append(b / m)
    ratios = np.array(ratios)
    print(f"{name:34s} applies to {n:5d} of {len(cubics):5d}; under-estimates {under}; bound/measured min {ratios.min():.4f} median {np.median(ratios):.3f} 90% {np.percentile(ratios, 90):.2f}" + (f"  WORST {worst}" if worst else ""))


rng = np.random.default_rng(5)
corpus = np.array(json.load(open(sys.argv[1])), dtype=float)
run("corpus (4000 sampled)", corpus[rng.choice(len(corpus), 4000, replace=False)])


def constructed(count, ratio, scale=100.0, far=False):
    out = []
    while len(out) < count:
        c0, q, c3 = rng.uniform(-scale, scale, (3, 2))
        if far:   # nearly parallel tangents: q far away
            mid = (c0 + c3) / 2
            q = mid + (q - mid) * 200
        l1, l2 = ratio()
        out.append(np.array([c0, c0 + l1 * (q - c0), c3 - l2 * (c3 - q), c3]))
    return out


run("handle ratios uniform in [0, 2]", constructed(3000, lambda: rng.uniform(0, 2, 2)))
run("ratios at the edges 0 and 2", constructed(2000, lambda: rng.choice([0.0, 2.0, 1e-9, 2 - 1e-9, 1.0], 2)))
run("ratios near 2/3", constructed(2000, lambda: 2 / 3 + rng.uniform(-0.05, 0.05, 2)))
run("ratios summing to 4/3", constructed(2000, lambda: (lambda x: (x, 4 / 3 - x))(rng.uniform(0, 4 / 3))))
run("tangent intersection far away", constructed(2000, lambda: rng.uniform(0, 0.02, 2), far=True))
run("scale 1e-3", constructed(1000, lambda: rng.uniform(0, 2, 2), scale=1e-3))
run("scale 1e6", constructed(1000, lambda: rng.uniform(0, 2, 2), scale=1e6))
zero = []
for _ in range(2000):
    c0, c2, c3 = rng.uniform(-100, 100, (3, 2))
    zero.append(np.array([c0, c0, c2, c3]) if rng.random() < 0.5 else np.array([c0, c2, c3, c3]))
run("zero-length handle", zero)
