#!/usr/bin/env python3
"""Check every mathematical statement of the note, numerically and, where
the statement is an identity of polynomials, exactly in rational arithmetic.

Usage: check_math.py

Each check prints what it checked and the largest deviation it found. Any
failure raises.
"""

import math
import random
from fractions import Fraction as F

import numpy as np

import cu2qu_math as M

rnd = random.Random(2026)
LINE = "{:<64s} {}"


def rc(scale=100.0):
    return complex(rnd.uniform(-scale, scale), rnd.uniform(-scale, scale))


def random_cubic(scale=100.0):
    return [rc(scale) for _ in range(4)]


def ok(name, value, bound):
    assert value <= bound, (name, value, bound)
    print(LINE.format(name, "max deviation %.2e" % value))


# --- 1. Structure of the split pieces ------------------------------------
worst_a = worst_b = worst_c = 0.0
for _ in range(400):
    c = random_cubic()
    d3 = M.third_difference(c)
    for n in (2, 3, 4, 5, 7, 12):
        d = d3 / n**3
        mids = []
        for k in range(n):
            p = M.split_cubic(c, k / n, (k + 1) / n)
            # every piece has the same third difference, d
            worst_a = max(worst_a, abs(M.third_difference(p) - d))
            q0, q1 = M.tangent_candidates(p)
            # the two candidates are d/2 apart
            worst_b = max(worst_b, abs((q0 - q1) - d / 2))
            mids.append((q0 + q1) / 2)
        for k in range(1, n):
            # halfway points of consecutive pieces straddle the split point
            worst_c = max(worst_c, abs((mids[k - 1] + mids[k]) / 2 - M.cubic_at(c, k / n)))
ok("each piece's third difference is d3 / n^3", worst_a, 1e-10)
ok("the two tangent candidates are d/2 apart", worst_b, 1e-10)
ok("midpoint of consecutive halfway points = split point", worst_c, 1e-10)

# --- 2. The error is d * phi_k(u), for any weights ------------------------
worst = 0.0
for _ in range(300):
    c = random_cubic()
    for n in (2, 3, 4, 6, 9):
        d = M.third_difference(c) / n**3
        for weights in (None, [rnd.uniform(-1.5, 1.5) for _ in range(n)]):
            quads = M.spline_quadratics(c, M.spline(c, n, weights))
            for k in range(n):
                for u in (0.0, 0.13, 0.5, 0.77, 1.0):
                    diff = M.quadratic_at(quads[k], u) - M.cubic_at(c, (k + u) / n)
                    worst = max(worst, abs(diff - d * M.phi(n, k, u, weights)))
ok("spline - cubic = d * phi_k(u), ramp and random weights", worst, 1e-10)

# --- 3. Closed forms for the ramp ----------------------------------------
worst_first = worst_mid = 0.0
for n in range(2, 60):
    g = (n - 2) / (4 * (n - 1))
    cc = (n - 3) / (4 * (n - 1))
    for i in range(101):
        u = i / 100
        worst_first = max(worst_first, abs(M.phi(n, 0, u) - (u * u * (1 + g) - u**3)))
        for k in range(1, n - 1):
            v = u - 0.5
            worst_mid = max(worst_mid, abs(M.phi(n, k, u) - ((n - 1 - 2 * k) / (4 * (n - 1)) + cc * v - v**3)))
ok("first piece: u^2 (1 + g) - u^3", worst_first, 1e-14)
ok("pieces in between: (n-1-2k)/(4(n-1)) + c v - v^3", worst_mid, 1e-14)

# Phi(n) against the largest |phi| found by dense evaluation
worst = 0.0
U = np.linspace(0.0, 1.0, 20001)
for n in list(range(2, 130)) + [200, 500, 1000, 2000]:
    pieces = range(n) if n < 130 else (0, 1, 2, n // 2, n - 3, n - 2, n - 1)
    best = max(float(np.abs(M.phi(n, k, U)).max()) for k in pieces)
    assert best <= M.error_factor(n) * (1 + 1e-12), (n, best, M.error_factor(n))
    worst = max(worst, M.error_factor(n) - best)
ok("Phi(n) = largest |phi_k(u)| (n = 2..129, 200..2000)", worst, 1e-8)

assert M.error_factor(2) == 4 / 27 and abs(M.error_factor(3) - 27 / 128) < 1e-16
assert M.limit(2) == 54.0 and M.limit(3) == 128.0
lim = 0.25 + math.sqrt(3) / 36
assert all(M.error_factor(n) < M.error_factor(n + 1) < lim for n in range(2, 5000))
assert all(M.limit(n) < M.limit(n + 1) for n in range(2, 20000))


def first_piece(n):
    g = (n - 2) / (4 * (n - 1))
    return (4 / 27) * (1 + g) * (1 + g) * (1 + g)


def second_piece(n):
    c = (n - 3) / (4 * (n - 1))
    w = c / 3
    return c + 2 * w * math.sqrt(w)


first_wins = [n for n in range(4, 400) if first_piece(n) > second_piece(n)]
assert first_wins == list(range(4, 55)), first_wins
print(LINE.format("Phi(2) = 4/27, Phi(3) = 27/128, increasing to 1/4 + sqrt(3)/36", "ok (%.6f)" % lim))
print(LINE.format("L(n) = n^3 / Phi(n) increases (n up to 20000)", "ok"))
print(LINE.format("first piece is the worst for n <= 54, the second from 55", "ok"))

# bounds on n
for _ in range(20000):
    ratio = 10 ** rnd.uniform(0, 7)  # |d3| / tolerance
    n = M.count(ratio, 1.0, 10**6) if ratio < 1e6 else 0
    if n:
        lo = (4 * ratio / 27) ** (1 / 3)
        hi = max(2, math.ceil((lim * ratio) ** (1 / 3)))
        assert lo <= n * (1 + 1e-12) and n <= hi, (ratio, n, lo, hi)
print(LINE.format("cbrt(4 r / 27) <= n <= ceil(cbrt(0.2981 r)),  r = |d3| / tol", "ok"))

# --- 4. One kept end tangent costs 4/27 -----------------------------------
# the two identities of the proof, exactly
for _ in range(200):
    n = rnd.randint(1, 9)
    s = [F(rnd.randint(-30, 30), rnd.randint(1, 12)) for _ in range(n)]

    def g(j):
        return F(0) if j in (0, n) else (s[j - 1] + s[j]) / 8

    def phi_exact(k, u):
        return g(k) * (1 - u) ** 2 + g(k + 1) * u * u + u * (1 - u) * (2 * u - 1 + s[k]) / 2

    for k in range(n):
        assert phi_exact(k, F(2, 3)) == F(4, 27) + 4 * g(k + 1) / 9 + g(k) / 9 + (s[k] - 1) / 9
        assert phi_exact(k, F(1, 3)) == -F(4, 27) + g(k + 1) / 9 + 4 * g(k) / 9 + (s[k] + 1) / 9
print(LINE.format("phi_k(2/3), phi_k(1/3) identities of the lower bound", "exact"))


def max_phi(n, weights, steps=400):
    u = np.linspace(0.0, 1.0, steps + 1)
    return max(float(np.abs(M.phi(n, k, u, weights)).max()) for k in range(n))


worst = 1.0
for _ in range(3000):
    n = rnd.randint(1, 8)
    w = [1.0] + [rnd.uniform(-3, 3) for _ in range(n - 1)]
    worst = min(worst, max_phi(n, w, 300))
assert worst >= 4 / 27 - 1e-9, worst
print(LINE.format("random weights with s_0 = 1: least max |phi| found", "%.6f >= 4/27 = %.6f" % (worst, 4 / 27)))
for n in (2, 4, 6, 8, 20):
    w = [1.0 if k % 2 == 0 else -1.0 for k in range(n)]
    assert abs(max_phi(n, w, 3000) - 4 / 27) < 1e-7
print(LINE.format("alternating weights, even n: max |phi| = 4/27", "ok"))
# n = 3 with both tangents kept: one free weight, best at 0
best = min((max_phi(3, [1.0, s1 / 200, -1.0], 600), s1 / 200) for s1 in range(-400, 401))
assert abs(best[1]) < 1e-12 and abs(best[0] - 27 / 128) < 1e-7, best
print(LINE.format("n = 3, both tangents kept: best middle weight 0, 27/128", "ok"))

# --- 5. The single quadratic ----------------------------------------------
# no quadratic with the cubic's end points is closer than sqrt(3)/36 |d3|
worst_ratio = 10.0
at_best = 0.0
for _ in range(400):
    c = random_cubic()
    d3 = abs(M.third_difference(c))
    q0, q1 = M.tangent_candidates(c)
    mid = (q0 + q1) / 2
    for trial in range(8):
        q = mid if trial == 0 else mid + rc(30.0)
        err = max(abs(M.quadratic_at((c[0], q, c[3]), i / 600) - M.cubic_at(c, i / 600)) for i in range(601))
        if trial == 0:
            at_best = max(at_best, abs(err / d3 - math.sqrt(3) / 36))
        worst_ratio = min(worst_ratio, err / d3)
assert worst_ratio >= math.sqrt(3) / 36 - 1e-6
ok("quadratic with the cubic's end points: error = sqrt(3)/36 |d3| at best", at_best, 1e-6)
print(LINE.format("... and never less, for any off-curve point", "least ratio %.6f, sqrt(3)/36 = %.6f" % (worst_ratio, math.sqrt(3) / 36)))

# the identity B(t) = Q(tau(t)) + rho(t) (a - z), exactly
for trial in range(4):
    c0, q, c3 = [(F(rnd.randint(-50, 50)), F(rnd.randint(-50, 50))) for _ in range(3)]
    a = (q[0] - c0[0], q[1] - c0[1])
    z = (c3[0] - q[0], c3[1] - q[1])
    for l1 in (F(0), F(1, 3), F(2, 3), F(5, 4), F(2)):
        for l2 in (F(0), F(1, 2), F(2, 3), F(7, 5), F(2)):
            c1 = (c0[0] + l1 * a[0], c0[1] + l1 * a[1])
            c2 = (c3[0] - l2 * z[0], c3[1] - l2 * z[1])
            p, r = l1 / 2 - F(1, 3), F(1, 3) - l2 / 2
            for i in range(9):  # degree 6 in t
                t = F(i, 8)
                e = (1 - t) * p + t * r
                tau = 3 * t * (1 - t) ** 2 * (l1 / 2) + 3 * t * t * (1 - t) * (1 - l2 / 2) + t**3
                assert tau == t + 3 * t * (1 - t) * e
                rho = 3 * t * t * (1 - t) ** 2 * (l1 + l2 - F(4, 3) + 3 * e * e)
                for x in (0, 1):
                    B = (1 - t) ** 3 * c0[x] + 3 * t * (1 - t) ** 2 * c1[x] + 3 * t * t * (1 - t) * c2[x] + t**3 * c3[x]
                    Q = (1 - tau) ** 2 * c0[x] + 2 * tau * (1 - tau) * q[x] + tau * tau * c3[x]
                    assert B - Q - rho * (a[x] - z[x]) == 0
print(LINE.format("B(t) = Q(tau(t)) + rho(t) (a - z)", "exact, on a grid past its degree"))

# kink-free interpolation: u x v of an interpolated joint
worst = 0.0
for _ in range(2000):
    u1, u2 = rc(), rc()
    r1, r2 = rnd.uniform(0.2, 5), rnd.uniform(0.2, 5)
    lam = rnd.uniform(0, 1)
    u = (1 - lam) * u1 + lam * u2
    v = (1 - lam) * r1 * u1 + lam * r2 * u2
    cross = u.real * v.imag - u.imag * v.real
    want = (1 - lam) * lam * (r2 - r1) * (u1.real * u2.imag - u1.imag * u2.real)
    worst = max(worst, abs(cross - want))
ok("joint of an instance: u x v = (1-L) L (r2 - r1) (u1 x u2)", worst, 1e-9)

# --- 6. More about the weights ---------------------------------------------
# no weights at all beat the mid-point spline: max |phi| >= sqrt(3)/36
V0 = math.sqrt(3) / 36
worst = 1.0
for _ in range(4000):
    n = rnd.randint(1, 8)
    w = [rnd.uniform(-2, 2) * rnd.choice((1, 0.1, 0.01)) for _ in range(n)]
    worst = min(worst, max_phi(n, w, 600))
assert worst >= V0 - 1e-6, worst
for n in (1, 2, 3, 5, 8, 13):
    assert abs(max_phi(n, [0.0] * n, 6000) - V0) < 1e-8
print(LINE.format("any weights: max |phi| >= sqrt(3)/36; equal for all-zero weights", "least found %.6f, sqrt(3)/36 = %.6f" % (worst, V0)))
# the two identities of that proof, at u = 1/2 -+ sqrt(3)/6
ua, ub = 0.5 - math.sqrt(3) / 6, 0.5 + math.sqrt(3) / 6
worst = 0.0
for _ in range(300):
    n = rnd.randint(1, 9)
    w = [rnd.uniform(-3, 3) for _ in range(n)]
    g = [0.0] + [(w[j - 1] + w[j]) / 8 for j in range(1, n)] + [0.0]
    for k in range(n):
        worst = max(worst, abs(M.phi(n, k, ub, w) - (g[k] * ua * ua + g[k + 1] * ub * ub + V0 + w[k] / 12)))
        worst = max(worst, abs(M.phi(n, k, ua, w) - (g[k] * ub * ub + g[k + 1] * ua * ua - V0 + w[k] / 12)))
ok("phi_k at u = 1/2 -+ sqrt(3)/6, identities of the sqrt(3)/36 bound", worst, 1e-12)


# five quadratics with both derivatives kept: best weights (1, s, 0, -s, -1)
def five(a):
    return max_phi(5, [1.0, a, 0.0, -a, -1.0], 20000)


lo, hi = -1.5, 0.5
for _ in range(80):
    m1 = lo + (hi - lo) * 0.381966
    m2 = lo + (hi - lo) * 0.618034
    if five(m1) < five(m2):
        hi = m2
    else:
        lo = m1
a5 = (lo + hi) / 2
assert abs(a5 - (-0.933394)) < 2e-5 and abs(five(a5) - M.best_error_factor(5)) < 2e-8, (a5, five(a5))
print(LINE.format("n = 5, both derivatives kept: best weights (1, s, 0, -s, -1)", "s = %.5f, max |phi| = %.6f" % (a5, five(a5))))
assert abs(max_phi(4, [1, 1 / 3, -1 / 3, -1], 20000) - M.error_factor(4)) < 1e-9
print(LINE.format("the ramp for n = 4 is (1, 1/3, -1/3, -1): %.4f against 4/27 = %.4f" % (M.error_factor(4), 4 / 27), "ok"))

# the best weights for every n up to 16, by linear programming (needs scipy)
try:
    import best_weights
except ImportError:
    print(LINE.format("best weights by linear programming", "skipped: no scipy"))
else:
    worst = 0.0
    for n in range(2, 17):
        factor, weights = best_weights.best(n)
        assert abs(max_phi(n, list(weights), 4000) - factor) < 1e-7, n
        worst = max(worst, abs(factor - M.best_error_factor(n)))
        # no random weights with both ends kept do better
        for _ in range(60):
            trial = [1.0] + [w + rnd.gauss(0, 0.03) for w in weights[1:-1]] + [-1.0]
            assert max_phi(n, trial, 400) >= factor - 1e-5, (n, trial)
    ok("best weights with both ends kept, n = 2..16: the table in cu2qu_math", worst, 1e-9)
    odd = [M.best_error_factor(n) for n in range(3, 16, 2)]
    assert all(x > y > 4 / 27 for x, y in zip(odd, odd[1:]))
    print(LINE.format("... for odd n: decreasing towards 4/27", ", ".join("%.4f" % x for x in odd)))

# keeping only the direction of the end tangent: the point may slide along it
cubic = (0j, 1 + 1j, 1 + 1j, 2 + 0j)
d3 = abs(M.third_difference(cubic))


def slid(c, n, mu):
    """The ramp spline with its first and last off-curve points moved along
    the cubic's end tangents by mu times the handles."""
    off = M.spline(c, n)
    off[0] += mu * (c[1] - c[0])
    off[-1] -= mu * (c[3] - c[2])
    return off


def equal_parameter_error(c, off):
    n = len(off)
    quads = M.spline_quadratics(c, off)
    return max(abs(M.quadratic_at(quads[k], i / 400) - M.cubic_at(c, (k + i / 400) / n)) for k in range(n) for i in range(401))


best = min(equal_parameter_error(cubic, slid(cubic, 2, mu / 1000)) for mu in range(-200, 201))
assert best < (4 / 27) * d3 / 8 * 0.95, best
single = equal_parameter_error(cubic, [1 + 1j])
assert abs(single - 0.25) < 1e-6 and single < (4 / 27) * d3
print(LINE.format("tangent kept, handle not: (0,0) (1,1) (1,1) (2,0) gets within", "%.5f < %.5f = (4/27)|d3|/8 with two; %.2f < %.3f with one" % (best, (4 / 27) * d3 / 8, single, (4 / 27) * d3)))
# ... but not when the handles are at right angles to the third difference
arch = (0j, 1j, 1 + 1j, 1 + 0j)
d3 = M.third_difference(arch)
assert abs((arch[1] - arch[0]).real * d3.real + (arch[1] - arch[0]).imag * d3.imag) < 1e-12
for n in (2, 3, 4):
    least = min(equal_parameter_error(arch, slid(arch, n, mu / 500)) for mu in range(-150, 151))
    assert least >= M.best_error_factor(n) * abs(d3) / n**3 * (1 - 1e-4), (n, least)  # sampled
print(LINE.format("... an arch with handles across d3: sliding never beats the best weights", "ok"))

# L grows: Phi grows more slowly than n^3
worst = 0.0
for n in range(2, 5000):
    first = lambda m: (4 / 27) * (1 + (m - 2) / (4 * (m - 1))) ** 3  # noqa: E731
    assert abs(first(n + 1) / first(n) - (1 + 1 / (n * (5 * n - 6))) ** 3) < 1e-12
    if n >= 4:
        second = lambda m: (m - 3) / (4 * (m - 1)) + 2 * ((m - 3) / (12 * (m - 1))) ** 1.5  # noqa: E731
        assert second(n + 1) / second(n) < (1 + 2 / (n * (n - 3))) ** 1.5 * (1 + 1e-12)
        assert (1 + 2 / (n * (n - 3))) ** 1.5 < (1 + 1 / n) ** 3
    assert M.error_factor(n + 1) / M.error_factor(n) < ((n + 1) / n) ** 3
print(LINE.format("Phi(n+1)/Phi(n) < ((n+1)/n)^3, and the two factors given for it", "ok, n up to 5000"))

# --- 7. The spline's points -------------------------------------------------
worst_on = worst_handle = worst_lin = 0.0
for _ in range(300):
    c = random_cubic()
    d3 = M.third_difference(c)
    for n in (2, 3, 4, 7):
        off = M.spline(c, n)
        quads = M.spline_quadratics(c, off)
        for k in range(1, n):
            g = (n - 2 * k) / (4 * (n - 1))
            worst_on = max(worst_on, abs(quads[k][0] - (M.cubic_at(c, k / n) + g * d3 / n**3)))
        worst_handle = max(worst_handle, abs((off[0] - c[0]) - 1.5 / n * (c[1] - c[0])))
        worst_handle = max(worst_handle, abs((c[3] - off[-1]) - 1.5 / n * (c[3] - c[2])))
        # converting commutes with interpolating
        c2 = random_cubic()
        lam = rnd.uniform(-0.5, 1.5)
        mixed = [(1 - lam) * a + lam * b for a, b in zip(c, c2)]
        want = [(1 - lam) * a + lam * b for a, b in zip(off, M.spline(c2, n))]
        worst_lin = max(worst_lin, max(abs(a - b) for a, b in zip(M.spline(mixed, n), want)))
ok("implied on-curve point k = C(k/n) + g_k d, g_k = (n-2k)/(4(n-1))", worst_on, 1e-10)
ok("end handles of the spline = 3/(2n) of the cubic's", worst_handle, 1e-10)
ok("spline of an interpolated cubic = interpolated splines", worst_lin, 1e-9)

# one quadratic: its handles are the cubic's divided by the fractions l1, l2
worst = 0.0
for _ in range(300):
    c = random_cubic()
    q = M.tangent_intersection(c)
    if q is None:
        continue
    l1, l2 = M.handle_fractions(c, q)
    worst = max(worst, abs((q - c[0]) * l1 - (c[1] - c[0])), abs((c[3] - q) * l2 - (c[3] - c[2])))
ok("one quadratic: (q - c0) l1 = c1 - c0 and (c3 - q) l2 = c3 - c2", worst, 1e-7)

# quarter circle of radius 1 drawn as a cubic: bound of sec. 7
kappa = 4 / 3 * (math.sqrt(2) - 1)
circle = (complex(1, 0), complex(1, kappa), complex(kappa, 1), complex(0, 1))
bound = M.single_quadratic_bound(circle)
parabola_mid = abs(M.quadratic_at((circle[0], complex(1, 1), circle[3]), 0.5)) - 1
assert abs(bound - 0.0607) < 5e-4 and abs(parabola_mid - (3 * math.sqrt(2) / 4 - 1)) < 1e-12
print(LINE.format("quarter circle: bound %.4f R; the parabola is %.4f R outside at its middle" % (bound, parabola_mid), "ok"))
# two masters whose curves are the same sets, and their instance halfway
m1c = (0j, 1 / 3 + 0j, 2 / 3 + 0j, 1 + 0j)
m1q = (0j, 0.1 + 0j, 1 + 0j)
m2c = (0j, 1j / 3, 2j / 3, 1j)
m2q = (0j, 0.5j, 1j)
for c, q in ((m1c, m1q), (m2c, m2q)):
    # both are the straight segment between the same end points
    for i in range(101):
        for z in (M.cubic_at(c, i / 100), M.quadratic_at(q, i / 100)):
            along = (z - c[0]) / (c[3] - c[0])
            assert abs(along.imag) < 1e-12 and -1e-12 <= along.real <= 1 + 1e-12
half_q = tuple((a + b) / 2 for a, b in zip(m1q, m2q))
assert half_q == (0j, 0.05 + 0.25j, 0.5 + 0.5j)
bulge = max(abs((M.quadratic_at(half_q, i / 1000) * (1 - 1j)).imag) / math.sqrt(2) for i in range(1001))
assert abs(bulge - 0.0707) < 1e-4
print(LINE.format("two masters, each quadratic on its cubic; halfway it is off by", "%.4f" % bulge))

# the offset of Theorem 4 vanishes at the middle when l1 + l2 = 4/3 - 3/16 (l1 - l2)^2
for _ in range(200):
    l1, l2 = rnd.uniform(0, 2), rnd.uniform(0, 2)
    e = 0.5 * (l1 / 2 - 1 / 3) + 0.5 * (1 / 3 - l2 / 2)
    assert abs((l1 + l2 - 4 / 3 + 3 * e * e) - (l1 + l2 - 4 / 3 + 3 / 16 * (l1 - l2) ** 2)) < 1e-12
print(LINE.format("rho(1/2) = 0 exactly when l1 + l2 = 4/3 - (3/16)(l1 - l2)^2", "ok"))


# the bound of Theorem 4 against the measured distance, for handle fractions
# all over its range, the ends 0 and 2 included
def set_distance(f, g):
    worst = 0.0
    for source, target in ((f, g), (g, f)):
        worst = max(worst, float(np.abs(source[::8, None] - target[None, :]).min(axis=1).max()))
    return worst


ts = np.linspace(0, 1, 2401)
least_slack = math.inf
largest_ratio = 0.0
fractions = [(0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (0.0, 4 / 3), (2 / 3, 2 / 3), (0.27, 0.97)]
fractions += [(rnd.uniform(0, 2), rnd.uniform(0, 2)) for _ in range(300)]
for l1, l2 in fractions:
    c0, q, c3 = rc(), rc(), rc()
    a, z = q - c0, c3 - q
    cubic = (c0, c0 + l1 * a, c3 - l2 * z, c3)
    bound = M.distance_bound(l1, l2, a, z)
    assert bound is not None
    cubic_points = np.array([M.cubic_at(cubic, t) for t in ts])
    quadratic_points = np.array([M.quadratic_at((c0, q, c3), t) for t in ts])
    measured = set_distance(cubic_points, quadratic_points)
    # the samples are within half a step of the curves
    step = max(np.abs(np.diff(cubic_points)).max(), np.abs(np.diff(quadratic_points)).max())
    assert measured <= bound + step, (l1, l2, measured, bound)
    least_slack = min(least_slack, bound - measured + step)
    if l1 == l2 == 2 / 3:
        assert bound < 1e-12
assert M.distance_bound(-0.01, 1, 1j, 1 + 0j) is None and M.distance_bound(1, 2.01, 1j, 1 + 0j) is None
print(LINE.format("the bound of Theorem 4 is never below the measured distance, l1, l2 in [0, 2]", "ok, 306 cubics"))

# a quadratic that passes at equal parameters has e1 and e2 within 2 sqrt(3)
# of the tolerance; so a tangent turned around takes a short handle
us = np.linspace(0, 1, 4001)
worst = 0.0
for _ in range(400):
    e1, e2 = rc(), rc()
    if rnd.random() < 0.1:
        e2 = -e1  # the case that reaches the bound
    largest = float(np.abs(3 * us * (1 - us) * ((1 - us) * e1 + us * e2)).max())
    worst = max(worst, max(abs(e1), abs(e2)) / largest)
assert 2 * math.sqrt(3) - 1e-3 < worst <= 2 * math.sqrt(3) * (1 + 1e-6), worst  # the largest error is sampled
turned = (0j, 2 + 0j, -2 + 0j, 8 - 2j)
q = M.tangent_intersection(turned)
assert q == -2 + 0j  # behind c0: the quadratic leaves c0 away from c1
e1 = (turned[0] + 2 * q) / 3 - turned[1]
e2 = (2 * q + turned[3]) / 3 - turned[2]
largest = float(np.abs(3 * us * (1 - us) * ((1 - us) * e1 + us * e2)).max())
assert largest <= 1 and abs(turned[1] - turned[0]) == 2
print(LINE.format("|e1|, |e2| <= 2 sqrt(3) x the largest equal-parameter error", "largest ratio %.4f; a handle of 2 turned around at error %.3f" % (worst, largest)))

# FontForge's observation, for the points where the end tangents of each
# piece meet: the joints are near the middle of them, not at it
lopsided = (0j, 10 + 120j, 300 + 150j, 320 + 0j)
worst_meeting = worst_halfway = 0.0
for n in (2, 3, 4, 8):
    pieces = [M.split_cubic(lopsided, k / n, (k + 1) / n) for k in range(n)]
    meet = [M.tangent_intersection(p) for p in pieces]
    half = [sum(M.tangent_candidates(p)) / 2 for p in pieces]
    for k in range(1, n):
        worst_meeting = max(worst_meeting, abs((meet[k - 1] + meet[k]) / 2 - pieces[k][0]))
        worst_halfway = max(worst_halfway, abs((half[k - 1] + half[k]) / 2 - pieces[k][0]))
assert worst_meeting > 3 and worst_halfway < 1e-10
print(LINE.format("joints are midway between halfway points; between tangent intersections only nearly", "off by up to %.2f units on an arch 320 wide" % worst_meeting))

# constants of the related work
assert abs(M.error_factor(2) / 8 - 1 / 54) < 1e-15  # Truong, Yuksel and Seiler's eq. 29
assert abs(1 / (12 * math.sqrt(3)) - math.sqrt(3) / 36) < 1e-15 and abs((12 * math.sqrt(3)) ** 2 - 432) < 1e-9
# Chebyshev: t^3 minus its best quadratic on [0, 1] is T3(2t - 1) / 32
cheb = 32 * us**3 - 48 * us**2 + 18 * us - 1
assert abs(np.abs(cheb).max() - 1) < 1e-12
assert [round(float(32 * t**3 - 48 * t**2 + 18 * t - 1)) for t in (0, 0.25, 0.75, 1)] == [-1, 1, -1, 1]
assert abs(math.sqrt(3) / 36 - 1 / (32 * math.cos(math.pi / 6) ** 3)) < 1e-15
print(LINE.format("1/54 = Phi(2)/8; 432 = (12 sqrt(3))^2; end points free: 1/32, kept: sqrt(3)/36", "ok"))
print("all checks passed")
