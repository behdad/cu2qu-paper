#!/usr/bin/env python3
"""Check that cu2qu-figures.js computes what cu2qu_math.py computes.

Usage: check_js.py cu2qu-figures.js

The figures are drawn by JavaScript, which has the formulas a second time.
This runs them under node and compares the values with the Python ones: the
error factor and the limits, the three counts of Figure 3, the error
functions, the points of a spline, the bound for one quadratic (from a cubic,
and from its handle fractions), and the error of the single quadratic at
equal parameters, by which Figure 1 says where one quadratic would do.
"""

import json
import math
import random
import subprocess
import sys

import cu2qu_math as M
import measure

NODE = r"""
const L = require(process.argv[1]);
const F = L.math;
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const out = {};
out.errorFactor = input.ns.map((n) => F.errorFactor(n));
out.limit = input.ns.map((n) => F.limit(n));
out.count = input.rs.map((r) => F.count(r, 1, 1000));
out.countBest = input.rs.map((r) => F.countBest(r));
out.countMidpoint = input.rs.map((r) => F.countMidpoint(r));
out.phi = input.phis.map(([n, k, u]) => F.phi(n, k, u));
out.splines = input.cubics.map(([c, n]) => F.spline(c, n));
out.bounds = input.cubics.map(([c]) => F.singleBound(c));
out.d3 = input.cubics.map(([c]) => F.thirdDifference(c));
out.at = input.cubics.map(([c]) => F.cubicAt(c, 0.3));
out.fractionBounds = input.fractions.map(([l1, l2, a, z]) => F.distanceBound(l1, l2, a, z));
out.singleErrors = input.cubics.map(([c]) => F.singleError(c));
// figure 1, as first drawn, and for a cubic that is a quadratic
const first = L.figures.construction.draw(L.figures.construction.state());
const flat = L.figures.construction.draw({ c: [[0, 0], [160, 200], [320, 200], [480, 0]], tolerance: 10 });
out.one = [[first.n, first.one, first.readout], [flat.n, flat.one, flat.readout]];
console.log(JSON.stringify(out));
"""


def n_midpoint(r):
    return max(1, math.ceil((r * r / 432.0) ** (1 / 6)))


def n_best(r):
    for n in range(2, 1000):
        if r <= n**3 / M.best_error_factor(n):
            return n
    return 0


def main():
    lib = sys.argv[1]
    rnd = random.Random(7)
    ns = list(range(2, 101)) + [150, 400, 1000]
    rs = [54, 54.0000001, 128, 128.0000001, 432, 20.78, 20.79] + [10 ** rnd.uniform(0.5, 5) for _ in range(3000)]
    phis = []
    for _ in range(2000):
        n = rnd.randint(2, 20)
        phis.append([n, rnd.randrange(n), rnd.random()])
    cubics = []
    for _ in range(500):
        c = [[rnd.uniform(-300, 300), rnd.uniform(-300, 300)] for _ in range(4)]
        cubics.append([c, rnd.randint(2, 9)])
    # cubics made from a quadratic, where the bound applies
    for _ in range(500):
        c0, q, c3 = [complex(rnd.uniform(-300, 300), rnd.uniform(-300, 300)) for _ in range(3)]
        l1, l2 = rnd.uniform(0, 2), rnd.uniform(0, 2)
        c1, c2 = c0 + l1 * (q - c0), c3 - l2 * (c3 - q)
        cubics.append([[[p.real, p.imag] for p in (c0, c1, c2, c3)], rnd.randint(2, 9)])
    fractions = [[0, 0], [0, 4 / 3], [2, 2], [2 / 3, 2 / 3], [0.27, 0.97], [2.0000001, 1], [1, -0.0000001]]
    fractions += [[rnd.uniform(-0.2, 2.2), rnd.uniform(-0.2, 2.2)] for _ in range(300)]
    fractions = [[l1, l2, [rnd.uniform(-300, 300), rnd.uniform(-300, 300)], [rnd.uniform(-300, 300), rnd.uniform(-300, 300)]] for l1, l2 in fractions]
    payload = json.dumps(dict(ns=ns, rs=rs, phis=phis, cubics=cubics, fractions=fractions))
    got = json.loads(subprocess.run(["node", "-e", NODE, lib], input=payload, capture_output=True, text=True, check=True).stdout)

    worst = 0.0
    for n, a, b in zip(ns, got["errorFactor"], got["limit"]):
        assert a == M.error_factor(n) and b == M.limit(n), (n, a, b)
    print("error factor and limit, n = 2..100, 150, 400, 1000: the same bits")
    for r, a, b, c in zip(rs, got["count"], got["countBest"], got["countMidpoint"]):
        assert (a, b, c) == (M.count(r, 1.0, 1000), n_best(r), n_midpoint(r)), (r, a, b, c)
    print("the three counts for %d values of |d3| / tolerance: the same" % len(rs))
    for (n, k, u), a in zip(phis, got["phi"]):
        worst = max(worst, abs(a - M.phi(n, k, u)))
    assert worst < 1e-15, worst
    print("phi at %d points: largest difference %.1e" % (len(phis), worst))
    worst = 0.0
    bounds = 0
    for (c, n), off, bound, d3, at in zip(cubics, got["splines"], got["bounds"], got["d3"], got["at"]):
        z = [complex(*p) for p in c]
        for a, b in zip(off, M.spline(z, n)):
            worst = max(worst, abs(complex(*a) - b))
        worst = max(worst, abs(complex(*d3) - M.third_difference(z)), abs(complex(*at) - M.cubic_at(z, 0.3)))
        want = M.single_quadratic_bound(z)
        assert (bound is None) == (want is None), (c, bound, want)
        if want is not None:
            bounds += 1
            assert abs(bound - want) <= 1e-9 * max(1.0, want), (c, bound, want)
    assert worst < 1e-9, worst
    print("spline points, third differences and points of %d cubics: largest difference %.1e" % (len(cubics), worst))
    print("the bound for one quadratic: the same for %d cubics, and not applicable for the same %d" % (bounds, len(cubics) - bounds))
    applies = 0
    for (l1, l2, a, z), bound in zip(fractions, got["fractionBounds"]):
        want = M.distance_bound(l1, l2, complex(*a), complex(*z))
        assert (bound is None) == (want is None), (l1, l2, bound, want)
        if want is not None:
            applies += 1
            assert abs(bound - want) <= 1e-9 * max(1.0, want), (l1, l2, bound, want)
    assert got["fractionBounds"][0] is not None and got["fractionBounds"][2] is not None and got["fractionBounds"][5] is None and got["fractionBounds"][6] is None
    print("the same bound from handle fractions: the same for %d pairs, the ends 0 and 2 among them, and none outside for the same %d" % (applies, len(fractions) - applies))
    worst = 0.0
    compared = 0
    for (c, n), error in zip(cubics, got["singleErrors"]):
        want = measure.single_equal_parameter_error([complex(*q) for q in c])
        assert (error is None) == (want is None), (c, error, want)
        if want:
            compared += 1
            worst = max(worst, abs(error - want) / want)
    # the figure samples the difference at 1001 parameters; measure.py solves for its largest value
    assert worst < 1e-4, worst
    print("the single quadratic's error at equal parameters, sampled, against its closed form: %d cubics, largest relative difference %.1e" % (compared, worst))
    (n1, one1, text1), (n2, one2, text2) = got["one"]
    assert (n1, one1) == (3, False) and "one would do" not in text1
    assert (n2, one2) == (2, True) and "2 quadratics; one would do here, which cu2qu tries first (section 5)." in text2
    print("figure 1 says where one quadratic would do, and not where it would not")
    print("cu2qu-figures.js agrees with cu2qu_math.py")


if __name__ == "__main__":
    main()
