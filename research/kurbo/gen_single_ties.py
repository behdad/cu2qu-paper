#!/usr/bin/env python3
"""Inputs where the SINGLE quadratic (tangent intersection) is at its limit:
the accuracy is the largest distance between that quadratic and the cubic,
give or take a few ulps. This is where the recursive fit test is at its
least stable.

Usage: gen_single_ties.py <out> [count]
Format as kdiff_py.py (S <tol> <8 coords>, floats as hex bits).
"""
import math, random, struct, sys

out = open(sys.argv[1], "w")
count = int(sys.argv[2]) if len(sys.argv) > 2 else 30000
rnd = random.Random(99)


def bits(x):
    return struct.pack(">d", float(x)).hex()


def max_error(e1, e2):
    """Largest |3u(1-u)^2 e1 + 3u^2(1-u) e2| over [0, 1]."""

    def f(u):
        a = 3 * u * (1 - u) * (1 - u)
        b = 3 * u * u * (1 - u)
        return math.hypot(a * e1[0] + b * e2[0], a * e1[1] + b * e2[1])

    best_u, best = 0.0, 0.0
    N = 400
    for i in range(1, N):
        v = f(i / N)
        if v > best:
            best, best_u = v, i / N
    lo, hi = best_u - 1 / N, best_u + 1 / N
    for _ in range(200):  # ternary search to full precision
        m1 = lo + (hi - lo) / 3
        m2 = hi - (hi - lo) / 3
        if f(m1) < f(m2):
            lo = m1
        else:
            hi = m2
    return f((lo + hi) / 2)


n = 0
while n < count:
    size = rnd.choice([40, 200, 1000])
    p0 = (rnd.randint(-size, size), rnd.randint(-size, size))
    q = (rnd.randint(-size, size), rnd.randint(-size, size))
    p3 = (rnd.randint(-size, size), rnd.randint(-size, size))
    cross = (q[0] - p0[0]) * (p3[1] - q[1]) - (q[1] - p0[1]) * (p3[0] - q[0])
    if cross == 0:
        continue
    # handles along the tangents, at fractions near 2/3, on a grid of 1/den
    den = rnd.choice([1, 2, 4, 8, 3, 6, 12])
    a = 2 / 3 + rnd.choice([-1, 1]) * rnd.choice([0.01, 0.03, 0.1, 0.2])
    b = a if rnd.random() < 0.5 else 2 / 3 + rnd.choice([-1, 1]) * rnd.choice([0.01, 0.03, 0.1, 0.2])
    if rnd.random() < 0.5:
        # snapped to a grid: the tangent intersection is then only near q
        c1 = (round((p0[0] + a * (q[0] - p0[0])) * den) / den, round((p0[1] + a * (q[1] - p0[1])) * den) / den)
        c2 = (round((p3[0] + b * (q[0] - p3[0])) * den) / den, round((p3[1] + b * (q[1] - p3[1])) * den) / den)
        # the intersection of the two handles, as a float
        ab = (c1[0] - p0[0], c1[1] - p0[1])
        cd = (p3[0] - c2[0], p3[1] - c2[1])
        pcd = ab[0] * cd[1] - ab[1] * cd[0]
        if pcd == 0:
            continue
        h = (ab[0] * (p0[1] - c2[1]) - ab[1] * (p0[0] - c2[0])) / pcd
        t = (c2[0] + cd[0] * h, c2[1] + cd[1] * h)
    else:
        c1 = (p0[0] + a * (q[0] - p0[0]), p0[1] + a * (q[1] - p0[1]))
        c2 = (p3[0] + b * (q[0] - p3[0]), p3[1] + b * (q[1] - p3[1]))
        t = q
    e1 = (p0[0] + (t[0] - p0[0]) * (2 / 3) - c1[0], p0[1] + (t[1] - p0[1]) * (2 / 3) - c1[1])
    e2 = (p3[0] + (t[0] - p3[0]) * (2 / 3) - c2[0], p3[1] + (t[1] - p3[1]) * (2 / 3) - c2[1])
    err = max_error(e1, e2)
    if not (err > 1e-6):
        continue
    for k in (-3, -2, -1, 0, 1, 2, 3):
        tol = err
        for _ in range(abs(k)):
            tol = math.nextafter(tol, math.inf if k > 0 else -math.inf)
        out.write("S %s %s\n" % (bits(tol), " ".join(bits(v) for p in (p0, c1, c2, p3) for v in p)))
    n += 1
out.close()
print("cubics:", n, "lines:", n * 7, file=sys.stderr)
