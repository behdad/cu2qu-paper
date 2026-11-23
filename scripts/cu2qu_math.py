"""The formulas of the cu2qu note, in one place.

Everything here is plain Python on complex numbers (a point is x + y*1j),
written to be read next to the note, not for speed.
"""

import math


# ---------------------------------------------------------------- the spline


def split_cubic(c, t0, t1):
    """The piece of the cubic c = (c0, c1, c2, c3) between parameters t0, t1,
    as four control points (blossoming / de Casteljau)."""

    def blossom(a, b, cc):
        c0, c1, c2, c3 = c
        p01 = c0 + (c1 - c0) * a
        p12 = c1 + (c2 - c1) * a
        p23 = c2 + (c3 - c2) * a
        p012 = p01 + (p12 - p01) * b
        p123 = p12 + (p23 - p12) * b
        return p012 + (p123 - p012) * cc

    return (
        blossom(t0, t0, t0),
        blossom(t0, t0, t1),
        blossom(t0, t1, t1),
        blossom(t1, t1, t1),
    )


def tangent_candidates(p):
    """The two candidates for the off-curve point of a quadratic for the
    cubic piece p: the one with the cubic's tangent at the start, and the one
    with its tangent at the end."""
    p0, p1, p2, p3 = p
    return p0 + 1.5 * (p1 - p0), p3 + 1.5 * (p2 - p3)


def spline(c, n, weights=None):
    """The spline of n >= 2 quadratics for the cubic c: the list of its n
    off-curve points. On-curve points are implied: c0, the midpoints of
    consecutive off-curve points, c3.

    weights[k] in [-1, 1] places the off-curve point of piece k between its
    two candidates: +1 is the start-tangent candidate, -1 the end-tangent
    one. cu2qu's ramp is 1 - 2k/(n-1).
    """
    if weights is None:
        weights = [1 - 2 * k / (n - 1) for k in range(n)]
    off = []
    for k in range(n):
        q0, q1 = tangent_candidates(split_cubic(c, k / n, (k + 1) / n))
        off.append((q0 + q1) / 2 + weights[k] * (q0 - q1) / 2)
    return off


def spline_quadratics(c, off):
    """The quadratics (start, control, end) of the spline with off-curve
    points `off` for the cubic c."""
    n = len(off)
    on = [c[0]] + [(off[k - 1] + off[k]) / 2 for k in range(1, n)] + [c[3]]
    return [(on[k], off[k], on[k + 1]) for k in range(n)]


def cubic_at(c, t):
    s = 1 - t
    return s * s * s * c[0] + 3 * s * s * t * c[1] + 3 * s * t * t * c[2] + t * t * t * c[3]


def quadratic_at(q, u):
    s = 1 - u
    return s * s * q[0] + 2 * s * u * q[1] + u * u * q[2]


def third_difference(c):
    return c[3] - 3 * c[2] + 3 * c[1] - c[0]


# ------------------------------------------------------- the error function


def phi(n, k, u, weights=None):
    """Quadratic k of the spline minus the cubic, at the parameter u of the
    piece, as a multiple of d = third_difference / n**3."""
    if weights is None:
        weights = [1 - 2 * j / (n - 1) for j in range(n)]

    def g(j):  # offset of on-curve point j from the cubic, in units of d
        return 0.0 if j in (0, n) else (weights[j - 1] + weights[j]) / 8

    s = weights[k]
    return g(k) * (1 - u) ** 2 + g(k + 1) * u * u + u * (1 - u) * (2 * u - 1 + s) / 2


def error_factor(n):
    """Phi(n): the largest |phi| over the spline of n >= 2 quadratics with
    cu2qu's ramp. The spline is |third difference| * Phi(n) / n**3 from the
    cubic at most, at equal parameters, and reaches it."""
    g = (n - 2) / (4 * (n - 1))
    error = (4 / 27) * (1 + g) * (1 + g) * (1 + g)  # first and last piece
    if n > 3:
        c = (n - 3) / (4 * (n - 1))
        w = c / 3
        error = max(error, c + 2 * w * math.sqrt(w))  # second and next to last
    return error


def limit(n):
    """L(n): the spline of n quadratics is within tolerance exactly when
    |third difference| <= tolerance * L(n)."""
    return n * n * n / error_factor(n)


def count(d3, tolerance, max_n=100):
    """The least n >= 2 whose spline is within tolerance, or 0."""
    for n in range(2, max_n + 1):
        if d3 <= tolerance * limit(n):
            return n
    return 0


# The best weights known that keep both end tangents (section 4 of the note):
# +1, -1, ... for even n; (1, 0, -1) for three; for odd n from five, one block
# (1, a, 0, -a, -1) with a = -0.9334 and pairs. check_math.py finds the block.
FIVE_BLOCK_FACTOR = 0.151879364835661
BEST_ODD_FACTOR = {3: 27 / 128, 5: FIVE_BLOCK_FACTOR}


def best_error_factor(n):
    if n % 2 == 0:
        return 4 / 27
    return 27 / 128 if n == 3 else FIVE_BLOCK_FACTOR


# --------------------------------------------------- the single quadratic


def tangent_intersection(c):
    """Where the two end tangents of the cubic meet, or None."""
    c0, c1, c2, c3 = c
    ab = c1 - c0
    cd = c3 - c2
    den = ab.real * cd.imag - ab.imag * cd.real
    if den == 0:
        return None
    h = (ab.real * (c0 - c2).imag - ab.imag * (c0 - c2).real) / den
    return c2 + cd * h


def handle_fractions(c, q):
    """l1, l2 with c1 = c0 + l1 (q - c0) and c2 = c3 - l2 (c3 - q)."""
    a = q - c[0]
    z = c[3] - q
    l1 = ((c[1] - c[0]) * a.conjugate()).real / abs(a) ** 2
    l2 = ((c[3] - c[2]) * z.conjugate()).real / abs(z) ** 2
    return l1, l2


def single_quadratic_bound(c):
    """An upper bound on the distance, as sets of points, between the cubic
    and the quadratic (c0, q, c3) of its end tangents; or None when the
    bound does not apply."""
    q = tangent_intersection(c)
    if q is None or q == c[0] or q == c[3]:
        return None
    l1, l2 = handle_fractions(c, q)
    if not (0 <= l1 <= 2 and 0 <= l2 <= 2):
        return None
    p = l1 / 2 - 1 / 3
    r = 1 / 3 - l2 / 2
    e2_max = max(p * p, r * r)
    e2_min = 0.0 if p * r <= 0 else min(p * p, r * r)
    g = l1 + l2 - 4 / 3
    m = max(abs(g + 3 * e2_max), abs(g + 3 * e2_min))
    a = q - c[0]
    z = c[3] - q
    return 3 / 16 * m * abs(a - z)
