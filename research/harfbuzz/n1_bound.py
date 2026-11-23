#!/usr/bin/env python3
"""How much does the simple geometric bound for a single quadratic buy?

For a cubic c0 c1 c2 c3 whose end tangents meet at q, with a = q - c0,
z = c3 - q, c1 = c0 + l1 a, c2 = c3 - l2 z:

    B(t) = Q(tau(t)) + rho(t) (a - z)
    tau(t)  = cubic Bezier function with control values 0, l1/2, 1 - l2/2, 1
    rho(t)  = 3 t^2 (1-t)^2 g2(t),   g2 = l1 + l2 - 4/3 + 3 e(t)^2
    e(t)    = (1-t) (l1/2 - 1/3) + t (1/3 - l2/2)

For 0 <= l1, l2 <= 2, tau maps [0,1] onto [0,1] inside [0,1], so every point
of the cubic is rho(t) (a - z) from a point of the quadratic's arc, and every
point of the arc is that far from a point of the cubic:

    Hausdorff distance <= (3/16) max|g2| |a - z|                (bound 1)

This script counts, on the corpus, how many cubics that accepts as a single
quadratic at each tolerance, against today's equal-parameter test and against
a measured distance, and checks that the bound is never below the measured
distance.

Usage: n1_bound.py cubics.json
"""
import json, math, sys
import numpy as np

cubics = np.array(json.load(open(sys.argv[1])), dtype=float)  # (N, 4, 2)
N = len(cubics)
c0, c1, c2, c3 = (cubics[:, i, :] for i in range(4))


def cross(u, v):
    return u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]


def dot(u, v):
    return u[:, 0] * v[:, 0] + u[:, 1] * v[:, 1]


# tangent intersection, as cu2qu's calc_intersect
ab = c1 - c0
cd = c3 - c2
den = cross(ab, cd)
with np.errstate(divide="ignore", invalid="ignore"):
    h = cross(ab, c0 - c2) / den
q = c2 + cd * h[:, None]
valid = np.isfinite(h) & (den != 0)

a = q - c0
z = c3 - q
with np.errstate(divide="ignore", invalid="ignore"):
    l1 = dot(c1 - c0, a) / dot(a, a)
    l2 = dot(c3 - c2, z) / dot(z, z)
ok = valid & (l1 >= 0) & (l1 <= 2) & (l2 >= 0) & (l2 <= 2)
p = l1 / 2 - 1 / 3
r = 1 / 3 - l2 / 2
e2max = np.maximum(p * p, r * r)
e2min = np.where(p * r <= 0, 0.0, np.minimum(p * p, r * r))
hi = 2 * (p - r) + 3 * e2max
lo = 2 * (p - r) + 3 * e2min
M = np.maximum(np.abs(hi), np.abs(lo))
amz = np.hypot((a - z)[:, 0], (a - z)[:, 1])
bound1 = np.where(ok, 3 / 16 * M * amz, np.inf)

# today's test: equal-parameter error of the tangent-intersection quadratic
# = max |3u(1-u)^2 e1 + 3u^2(1-u) e2|, by dense sampling (measurement only)
e1 = c0 + (q - c0) * (2 / 3) - c1
e2 = c3 + (q - c3) * (2 / 3) - c2
u = np.linspace(0, 1, 2001)
w1 = 3 * u * (1 - u) ** 2
w2 = 3 * u * u * (1 - u)
ep = np.full(N, np.inf)
idx = np.where(valid)[0]
for s in range(0, len(idx), 4000):
    i = idx[s : s + 4000]
    ex = e1[i, 0][:, None] * w1 + e2[i, 0][:, None] * w2
    ey = e1[i, 1][:, None] * w1 + e2[i, 1][:, None] * w2
    ep[i] = np.sqrt(ex * ex + ey * ey).max(axis=1)

# measured two-sided distance between the cubic and the quadratic, by dense
# sampling of both (800 points each) and nearest sample; good to about 1e-3 of
# the curve's size, which is enough to count and to catch a wrong bound
def measured(i):
    t = np.linspace(0, 1, 801)
    B = ((1 - t) ** 3)[:, None] * c0[i] + (3 * t * (1 - t) ** 2)[:, None] * c1[i] + (3 * t * t * (1 - t))[:, None] * c2[i] + (t**3)[:, None] * c3[i]
    Q = ((1 - t) ** 2)[:, None] * c0[i] + (2 * t * (1 - t))[:, None] * q[i] + (t * t)[:, None] * c3[i]
    d = np.sqrt(((B[:, None, :] - Q[None, :, :]) ** 2).sum(axis=2))
    return max(d.min(axis=1).max(), d.min(axis=0).max())


print(f"cubics {N}; tangents meet: {valid.sum()}; handle ratios in [0, 2]: {ok.sum()}")
for tol in (1.0, 0.5, 0.1):
    acc_ep = (ep <= tol).sum()
    acc_b1 = (bound1 <= tol).sum()
    acc_either = ((ep <= tol) | (bound1 <= tol)).sum()
    print(f"tol {tol}: accepted as one quadratic: equal-parameter test {acc_ep} ({acc_ep / N:.1%}); bound 1 {acc_b1} ({acc_b1 / N:.1%}); either {acc_either} ({acc_either / N:.1%})")

# safety and tightness of bound 1 on a sample of the cubics it applies to
rng = np.random.default_rng(1)
sample = rng.choice(np.where(ok & np.isfinite(bound1) & (bound1 < 50))[0], size=6000, replace=False)
ratios = []
under = 0
for i in sample:
    m = measured(i)
    size = max(np.abs(cubics[i]).max(), 1.0)
    if bound1[i] < m - 2e-3 * size:      # sampling resolution of the instrument
        under += 1
    if m > 1e-9:
        ratios.append(bound1[i] / m)
ratios = np.array(ratios)
print(f"bound 1 against measured distance on {len(sample)} cubics: under-estimates {under}; ratio bound/measured: min {ratios.min():.4f}, 25% {np.percentile(ratios, 25):.3f}, median {np.median(ratios):.3f}, 75% {np.percentile(ratios, 75):.3f}, 95% {np.percentile(ratios, 95):.2f}")

# how many of the cubics truly within tolerance of the tangent-intersection
# quadratic does bound 1 catch? (measured on a sample)
for tol in (1.0, 0.5):
    sample2 = rng.choice(np.where(valid)[0], size=4000, replace=False)
    within = caught = caught_ep = 0
    for i in sample2:
        m = measured(i)
        if m <= tol:
            within += 1
            caught += bound1[i] <= tol
            caught_ep += ep[i] <= tol
    print(f"tol {tol}: of {len(sample2)} sampled cubics, {within} are within tolerance of their tangent-intersection quadratic by measurement; bound 1 accepts {caught} of those, the equal-parameter test {caught_ep}")
