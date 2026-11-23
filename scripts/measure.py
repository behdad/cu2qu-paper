#!/usr/bin/env python3
"""Measure the conversion schemes of the note on a corpus of font cubics.

Usage: measure.py data/cubics.json [--hausdorff] > data/measure.txt

Everything here is computed from the formulas in cu2qu_math.py, independent
of fontTools, kurbo and HarfBuzz (compare_revs.py checks fontTools against
it). With --hausdorff it also measures the distance between each cubic and
its splines as sets of points, which takes a few minutes.
"""

import json
import math
import sys

import numpy as np

import cu2qu_math as M

TOLERANCES = (1.0, 0.5, 0.1)
SQRT3_36 = math.sqrt(3) / 36


def load(path):
    return [tuple(complex(x, y) for x, y in c) for c in json.load(open(path))]


def dot(a, b):
    return a.real * b.real + a.imag * b.imag


def cross(a, b):
    return a.real * b.imag - a.imag * b.real


# ---------------------------------------------------------------- one quadratic


def ft_intersection(c):
    """Where the end tangents meet, as fontTools' calc_intersect finds it."""
    c0, c1, c2, c3 = c
    ab = c1 - c0
    cd = c3 - c2
    p = ab * 1j
    den = dot(p, cd)
    if den == 0:
        if c1 == c2 and (c0 == c1 or c2 == c3):
            return c1
        return None
    h = dot(p, c0 - c2) / den
    return c2 + cd * h


def single_equal_parameter_error(c):
    """The largest distance, at equal parameters, between the cubic and the
    quadratic (c0, q, c3) of its end tangents; None if there is none.

    The difference is 3u(1-u)^2 e1 + 3u^2(1-u) e2; its length is largest at a
    root of 3c u^3 + (5b - 2c) u^2 + (2a - 3b) u - a, with a = |e1|^2,
    b = e1.(e2 - e1), c = |e2 - e1|^2.
    """
    q = ft_intersection(c)
    if q is None:
        return None
    e1 = (c[0] + 2 * q) / 3 - c[1]
    e2 = (2 * q + c[3]) / 3 - c[2]
    a = dot(e1, e1)
    b = dot(e1, e2 - e1)
    cc = dot(e2 - e1, e2 - e1)
    if a == 0 and cc == 0:
        return 0.0
    best = 0.0
    for u in np.roots([3 * cc, 5 * b - 2 * cc, 2 * a - 3 * b, -a]):
        if abs(u.imag) < 1e-9 and 0 < u.real < 1:
            u = u.real
            best = max(best, abs(3 * u * (1 - u) * ((1 - u) * e1 + u * e2)))
    return best


def hb_single_bound(c):
    """HarfBuzz's bound on the distance, as sets of points, between the cubic
    and the quadratic of its end tangents (hb_gpu_cubic_approx_quadratic);
    None where it gives up. Returns (bound, q)."""
    c0, c1, c2, c3 = c
    a_dir = c1 - c0
    d_dir = c3 - c2
    if a_dir == 0:
        a_dir = c2 - c0
    if d_dir == 0:
        d_dir = c3 - c1
    p = complex(-a_dir.imag, a_dir.real)
    den = dot(p, d_dir)
    if abs(den) < 1e-12:
        return None
    h = dot(p, c0 - c2) / den
    q = c2 + d_dir * h
    a = q - c0
    z = c3 - q
    aa = dot(a, a)
    zz = dot(z, z)
    if not (aa > 0 and zz > 0):
        return None
    l1 = dot(c1 - c0, a) / aa
    l2 = dot(c3 - c2, z) / zz
    if not (0 <= l1 <= 2 and 0 <= l2 <= 2):
        return None
    pp = l1 * 0.5 - 1 / 3
    r = 1 / 3 - l2 * 0.5
    e2_max = max(pp * pp, r * r)
    e2_min = 0.0 if pp * r <= 0 else min(pp * pp, r * r)
    g = l1 + l2 - 4 / 3
    m = max(abs(g + 3 * e2_max), abs(g + 3 * e2_min))
    return 3 / 16 * m * abs(a - z), q


def hb_draws_line(c, tol):
    """HarfBuzz's GPU library draws a cubic as a line, before its converter is
    reached, when both control points are within the tolerance of the line
    through the end points (hb_gpu_draw_t::acc_cubic_to)."""
    c0, c1, c2, c3 = c
    chord = c3 - c0
    length_sq = dot(chord, chord)
    if not length_sq > 0:
        return False
    cross1 = cross(c1 - c0, chord)
    cross2 = cross(c2 - c0, chord)
    return cross1 * cross1 / length_sq <= tol * tol and cross2 * cross2 / length_sq <= tol * tol


# ------------------------------------------------------------ numbers of pieces


def n_ramp(d3, tol):
    return M.count(d3, tol, 1000)


def n_midpoint(d3, tol):
    """kurbo's CubicBez::to_quads and lyon: the least n with
    sqrt(3)/36 |d3| / n^3 <= tol, computed as kurbo does."""
    return max(1, math.ceil((d3 * d3 / (432.0 * tol * tol)) ** (1 / 6)))


def n_best_weights(d3, tol):
    """The least n >= 2 for the best weights that keep both end handles."""
    for n in range(2, 1000):
        if d3 <= tol * n * n * n / M.best_error_factor(n):
            return n
    return 0


# ----------------------------------------------------------------- end tangents


def angle(a, b):
    """Angle between two directions, in degrees; None if either is zero."""
    if a == 0 or b == 0:
        return None
    return abs(math.degrees(math.atan2(cross(a, b), dot(a, b))))


def end_directions(c):
    """The directions in which the cubic leaves its start and arrives at its
    end (taken from the next control point when a handle has no length)."""
    c0, c1, c2, c3 = c
    start = c1 - c0 or c2 - c0 or c3 - c0
    end = c3 - c2 or c3 - c1 or c3 - c0
    return start, end


def midpoint_end_directions(c, n):
    """The same for the mid-point spline of n quadratics."""
    d3 = M.third_difference(c)
    start = 1.5 * (c[1] - c[0]) / n - d3 / (4 * n**3)
    end = 1.5 * (c[3] - c[2]) / n - d3 / (4 * n**3)
    return start, end


# ---------------------------------------------------- distance as sets of points


def power_cubic(c):
    c0, c1, c2, c3 = c
    return np.array([[c0, 3 * (c1 - c0), 3 * (c0 - 2 * c1 + c2), c3 - 3 * c2 + 3 * c1 - c0]])


def power_quadratics(quads):
    return np.array([[a, 2 * (o - a), a - 2 * o + b, 0] for a, o, b in quads])


def evaluate(coef, u):
    return coef[..., 0] + u * (coef[..., 1] + u * (coef[..., 2] + u * coef[..., 3]))


def derivative(coef, u):
    return coef[..., 1] + u * (2 * coef[..., 2] + u * 3 * coef[..., 3])


def second_derivative(coef, u):
    return 2 * coef[..., 2] + u * 6 * coef[..., 3]


def nearest_distance(points, coef, samples):
    """Distance from each point to the curve made of the polynomial pieces
    coef (pieces x 4, power basis, parameter 0..1 each): the best of a grid of
    samples, refined by Newton's method in that piece and in its neighbours."""
    pieces = coef.shape[0]
    grid = np.linspace(0.0, 1.0, samples + 1)
    curve = evaluate(coef[:, None, :], grid[None, :])  # pieces x samples+1
    d2 = np.abs(points[:, None] - curve.reshape(1, -1)) ** 2
    flat = d2.argmin(axis=1)
    piece = flat // (samples + 1)
    u = grid[flat % (samples + 1)]
    best = np.sqrt(d2[np.arange(len(points)), flat])
    for shift, start in ((0, None), (-1, 1.0), (1, 0.0)):
        i = piece + shift
        ok = (i >= 0) & (i < pieces)
        if not ok.any():
            continue
        i = np.where(ok, i, piece)
        k = coef[i]
        v = u.copy() if start is None else np.full(len(points), start)
        for _ in range(8):
            delta = evaluate(k, v) - points
            first = derivative(k, v)
            g = (delta.conj() * first).real
            gg = np.abs(first) ** 2 + (delta.conj() * second_derivative(k, v)).real
            step = np.where(gg > 0, g / np.where(gg > 0, gg, 1.0), 0.0)
            v = np.clip(v - step, 0.0, 1.0)
        dist = np.abs(evaluate(k, v) - points)
        best = np.where(ok, np.minimum(best, dist), best)
    return best


def directed_distance(source, target, per_piece=32):
    """The largest distance from a point of the source curve to the target
    curve, both given as pieces in the power basis."""
    pieces = source.shape[0]
    u = np.linspace(0.0, 1.0, per_piece * max(pieces, target.shape[0]) // pieces + 1)
    points = evaluate(source[:, None, :], u[None, :]).reshape(-1)
    d = nearest_distance(points, target, 24 * max(1, pieces // target.shape[0]))
    i = int(d.argmax())
    top = float(d[i])
    # parabola through the largest sample and its neighbours on the same piece
    step = len(u)
    if i % step not in (0, step - 1):
        a, b, c = d[i - 1], d[i], d[i + 1]
        den = a - 2 * b + c
        if den < 0:
            top = float(b - (a - c) ** 2 / (8 * den))
    return top


def hausdorff(c, quads):
    cubic = power_cubic(c)
    spline = power_quadratics(quads)
    return max(directed_distance(cubic, spline), directed_distance(spline, cubic))


def ramp_quadratics(c, n):
    """The quadratics cu2qu makes for n: the single one of the end tangents,
    or the spline of n. None if there is no single quadratic."""
    if n == 1:
        q = ft_intersection(c)
        if q is None or not (math.isfinite(q.real) and math.isfinite(q.imag)):
            return None
        return [(c[0], q, c[3])]
    return M.spline_quadratics(c, M.spline(c, n))


# ------------------------------------------------------------------------ report


def mean(values):
    return sum(values) / len(values)


def percentile(values, p):
    values = sorted(values)
    return values[min(len(values) - 1, int(p / 100 * len(values)))]


def main():
    path = sys.argv[1]
    with_hausdorff = "--hausdorff" in sys.argv
    cubics = load(path)
    total = len(cubics)
    out = {"cubics": total}
    print("corpus: %d cubics (%s)" % (total, path))

    d3 = [abs(M.third_difference(c)) for c in cubics]
    print(
        "|d3| percentiles 50 / 90 / 99 / max: %.1f / %.1f / %.1f / %.1f"
        % (percentile(d3, 50), percentile(d3, 90), percentile(d3, 99), max(d3))
    )
    zero = sum(1 for c in cubics if c[0] == c[1] or c[2] == c[3])
    print("a control point on its end point: %d (%.2f%%)" % (zero, 100 * zero / total))

    single = [single_equal_parameter_error(c) for c in cubics]
    hb = [hb_single_bound(c) for c in cubics]

    print("\none quadratic")
    print("%-58s %10s %10s %10s" % ("tolerance", *("%.1f" % t for t in TOLERANCES)))
    rows = [
        ("fits at equal parameters (cu2qu's test)", lambda e, b, t: e is not None and e <= t),
        ("ruled out by |d3| > 36/sqrt(3) x tolerance", None),
        ("within tolerance by the bound on the distance of sec. 7", lambda e, b, t: b is not None and b[0] <= t),
    ]
    for name, test in rows:
        if test is None:
            values = [sum(1 for x in d3 if x * SQRT3_36 > t) for t in TOLERANCES]
        else:
            values = [sum(1 for e, b in zip(single, hb) if test(e, b, t)) for t in TOLERANCES]
        print("%-58s %10s %10s %10s" % (name, *("%d (%.1f%%)" % (v, 100 * v / total) for v in values)))
        out[name] = values
    # the bound can never rule out a quadratic that fits
    bad = sum(1 for e, x in zip(single, d3) for t in TOLERANCES if e is not None and e <= t and x * SQRT3_36 > t * (1 + 1e-12))
    print("fits although ruled out: %d" % bad)
    assert bad == 0

    print("\nquadratics per cubic")
    print("%-58s %10s %10s %10s" % ("tolerance", *("%.1f" % t for t in TOLERANCES)))
    counts = {}

    def scheme(name, fn):
        counts[name] = [[fn(i, t) for i in range(total)] for t in TOLERANCES]
        print("%-58s %10.6f %10.6f %10.6f" % (name, *(mean(v) for v in counts[name])))
        out[name] = [mean(v) for v in counts[name]]

    def is_single(i, t):
        return single[i] is not None and single[i] <= t

    scheme("cu2qu: one if it fits at equal parameters, else the ramp", lambda i, t: 1 if is_single(i, t) else n_ramp(d3[i], t))
    scheme("the same with the best weights (sec. 4)", lambda i, t: 1 if is_single(i, t) else n_best_weights(d3[i], t))
    scheme("HarfBuzz: one by the bound of sec. 7, else the ramp", lambda i, t: 1 if hb[i] is not None and hb[i][0] <= t else n_ramp(d3[i], t))
    scheme("mid-point (kurbo to_quads, lyon)", lambda i, t: n_midpoint(d3[i], t))

    # HarfBuzz itself never gives its converter a cubic that is nearly straight
    lines = [[hb_draws_line(c, t) for c in cubics] for t in TOLERANCES]
    print("\nHarfBuzz draws a line where both control points are within tolerance of the chord's line")
    print("%-58s %10s %10s %10s" % ("cubics", *("%d (%.1f%%)" % (sum(v), 100 * sum(v) / total) for v in lines)))
    hb_name = "HarfBuzz: one by the bound of sec. 7, else the ramp"
    with_lines = [mean([1 if line else n for line, n in zip(v, w)]) for v, w in zip(lines, counts[hb_name])]
    print("%-58s %10.6f %10.6f %10.6f" % ("segments per cubic with those as one line", *with_lines))
    out["hb lines"] = [sum(v) for v in lines]

    name = "cu2qu: one if it fits at equal parameters, else the ramp"
    print("\ncu2qu, share of cubics by number of quadratics (1 / 2 / 3 / 4 / 5 and more), percent")
    for t, v in zip(TOLERANCES, counts[name]):
        share = [100 * sum(1 for n in v if n == k) / total for k in (1, 2, 3, 4)]
        share.append(100 - sum(share))
        print("  tolerance %.1f: %s   (most: %d)" % (t, " / ".join("%.1f" % s for s in share), max(v)))
        out["share %.1f" % t] = share
    for t, v in zip(TOLERANCES, counts["mid-point (kurbo to_quads, lyon)"]):
        share = [100 * sum(1 for n in v if n == k) / total for k in (1, 2, 3, 4)]
        share.append(100 - sum(share))
        print("  mid-point %.1f: %s   (most: %d)" % (t, " / ".join("%.1f" % s for s in share), max(v)))

    print("\nwho can gain from better weights: cubics with |d3| > 128 x tolerance")
    for t in TOLERANCES:
        v = sum(1 for x in d3 if x > 128 * t)
        print("  tolerance %.1f: %d (%.1f%%)" % (t, v, 100 * v / total))

    print("\nmid-point spline: angle between the cubic's end tangents and the spline's")
    for t, v in zip(TOLERANCES, counts["mid-point (kurbo to_quads, lyon)"]):
        worst = []
        for c, n in zip(cubics, v):
            want = end_directions(c)
            got = midpoint_end_directions(c, n)
            angles = [angle(w, g) for w, g in zip(want, got)]
            angles = [a for a in angles if a is not None]
            worst.append(max(angles) if angles else 0.0)
        over = [sum(1 for a in worst if a > lim) for lim in (0.1, 1, 5)]
        print(
            "  tolerance %.1f: more than 0.1 / 1 / 5 degrees on %s cubics; median %.2f, mean %.2f degrees"
            % (t, " / ".join("%d (%.1f%%)" % (o, 100 * o / total) for o in over), percentile(worst, 50), mean(worst))
        )
        out["midpoint angle %.1f" % t] = [100 * o / total for o in over]

    if with_hausdorff:
        measure_hausdorff(cubics, d3, single, hb, counts[name], out)

    with open(path.replace(".json", ".measure.json"), "w") as f:
        json.dump(out, f, indent=1)


def measure_hausdorff(cubics, d3, single, hb, eq_counts, out):
    """For every cubic, the distance as sets of points to the single quadratic
    and to the ramp spline of every n up to what cu2qu uses at the finest
    tolerance; and to HarfBuzz's single quadratic, which differs from cu2qu's
    where a control point is on its end point."""
    total = len(cubics)
    finest = eq_counts[-1]
    table = []
    hb_distance = []
    for i, c in enumerate(cubics):
        row = {}
        for n in range(1, max(finest[i], 2) + 1):
            quads = ramp_quadratics(c, n)
            row[n] = math.inf if quads is None else hausdorff(c, quads)
        table.append(row)
        hb_distance.append(None if hb[i] is None else hausdorff(c, [(c[0], hb[i][1], c[3])]))
        if i % 10000 == 0:
            print("  ... %d" % i, file=sys.stderr)

    print("\ndistance as sets of points (Hausdorff), measured")
    print("%-58s %10s %10s %10s" % ("tolerance", *("%.1f" % t for t in TOLERANCES)))
    least = []
    for t in TOLERANCES:
        v = []
        for i, row in enumerate(table):
            n = next((n for n in sorted(row) if row[n] <= t), None)
            assert n is not None, (i, t, row)
            v.append(n)
        least.append(v)
    print("%-58s %10.6f %10.6f %10.6f" % ("cu2qu's constructions, least n within tolerance", *(mean(v) for v in least)))
    out["least n by distance"] = [mean(v) for v in least]
    fill = []
    for t, v in zip(TOLERANCES, eq_counts):
        fill.append(mean([table[i][n] / t for i, n in enumerate(v)]))
    print("%-58s %10.6f %10.6f %10.6f" % ("cu2qu's output: distance / tolerance, mean", *fill))
    out["fill"] = fill
    one = [sum(1 for row in table if row[1] <= t) for t in TOLERANCES]
    print("%-58s %10s %10s %10s" % ("one quadratic within tolerance", *("%d (%.1f%%)" % (v, 100 * v / total) for v in one)))
    out["one by distance"] = one
    # how much of what the equal-parameter test costs is in cubics that one
    # quadratic would do for
    share = []
    for v, w in zip(eq_counts, least):
        saved = sum(a - b for a, b in zip(v, w))
        share.append(100 * sum(a - b for a, b in zip(v, w) if b == 1) / saved)
    print("%-58s %10.1f %10.1f %10.1f" % ("share of that saving in cubics that take one, percent", *share))

    # HarfBuzz's own single quadratic
    one_hb = [sum(1 for h in hb_distance if h is not None and h <= t) for t in TOLERANCES]
    print("%-58s %10s %10s %10s" % ("HarfBuzz's quadratic within tolerance", *("%d (%.1f%%)" % (v, 100 * v / total) for v in one_hb)))
    # HarfBuzz as it is for one quadratic, with its splines accepted by distance
    beyond = []
    for t, v in zip(TOLERANCES, eq_counts):
        w = []
        for i, row in enumerate(table):
            if hb[i] is not None and hb[i][0] <= t:
                w.append(1)
            else:
                w.append(next((n for n in sorted(row) if n >= 2 and row[n] <= t), v[i]))
        beyond.append(mean(w))
    print("%-58s %10.6f %10.6f %10.6f" % ("HarfBuzz, its splines accepted by distance", *beyond))

    # the bound of sec. 7 must never be below the measured distance
    worst = 0.0
    ratios = []
    for b, measured in zip(hb, hb_distance):
        if b is None:
            continue
        ratios.append(b[0] / measured if measured > 1e-9 else 1.0)
        worst = max(worst, measured - b[0])
    print("bound of sec. 7 against the measured distance: below it by at most %.2e units; bound / measured median %.2f, 90th percentile %.2f (%d cubics)" % (worst, percentile(ratios, 50), percentile(ratios, 90), len(ratios)))
    # at equal parameters the distance is an upper bound of the distance as sets
    over = 0.0
    for i, row in enumerate(table):
        for n in row:
            if n >= 2:
                eq = d3[i] * M.error_factor(n) / n**3
                over = max(over, row[n] - eq)
    print("measured distance above the equal-parameter distance by at most %.2e units (n >= 2)" % over)


def hausdorff_single(row):
    return row[1] if math.isfinite(row[1]) else None


if __name__ == "__main__":
    main()
