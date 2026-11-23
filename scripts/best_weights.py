#!/usr/bin/env python3
"""The best weights for a spline that keeps both end handles.

Usage: best_weights.py [MAX_N]

For each n from 2, find the weights s_0 = 1, s_1, ..., s_{n-1} = -1 that
make the largest |phi_k(u)| over the spline least. phi is linear in the
weights, so this is a linear program: minimise M subject to
-M <= phi_k(u) <= M at chosen parameters u. The parameters start as a grid;
after each solution the exact largest |phi_k| of every piece is found (each
phi_k is a cubic in u, so at an end or at a root of its derivative), and the
parameters where it exceeds M are added, until none does.

Prints the least largest value (the error factor), the limit n^3 / factor,
and the weights; and, with --table, the dictionary that cu2qu_math.py holds.
The factor printed is the exact largest |phi| of the weights found, so those
weights do reach it; that no weights do better is as good as the solver's
tolerance, about nine digits, which is why this stops at sixteen: beyond that
the best factor for an odd n is within 1e-7 of 4/27. Needs scipy.
"""

import math
import sys

import numpy as np
from scipy.optimize import linprog


def row(n, k, u):
    """phi_k(u) as constant + coefficients . (s_0, ..., s_{n-1})."""
    c = np.zeros(n)
    if k >= 1:  # g_k = (s_{k-1} + s_k) / 8
        c[k - 1] += (1 - u) ** 2 / 8
        c[k] += (1 - u) ** 2 / 8
    if k + 1 <= n - 1:  # g_{k+1} = (s_k + s_{k+1}) / 8
        c[k] += u * u / 8
        c[k + 1] += u * u / 8
    c[k] += u * (1 - u) / 2
    return u * (1 - u) * (2 * u - 1) / 2, c


def extremes(n, k, s):
    """The parameters where |phi_k| can be largest: the ends, and the roots
    of the derivative of the cubic phi_k."""
    g0 = 0.0 if k == 0 else (s[k - 1] + s[k]) / 8
    g1 = 0.0 if k == n - 1 else (s[k] + s[k + 1]) / 8
    # phi = g0 (1-u)^2 + g1 u^2 + (u - u^2)(2u - 1 + s_k) / 2
    #     = g0 + (-2 g0 + (s_k - 1)/2) u + (g0 + g1 + (3 - s_k)/2) u^2 - u^3
    b, c = -2 * g0 + (s[k] - 1) / 2, g0 + g1 + (3 - s[k]) / 2
    us = [0.0, 1.0]
    disc = 4 * c * c + 12 * b  # of -3 u^2 + 2 c u + b
    if disc >= 0:
        for sign in (-1, 1):
            u = (2 * c + sign * math.sqrt(disc)) / 6
            if 0 < u < 1:
                us.append(u)
    return us


def value(n, k, u, s):
    const, coef = row(n, k, u)
    return const + float(coef @ s)


def largest(n, s):
    return max(abs(value(n, k, u, s)) for k in range(n) for u in extremes(n, k, s))


def best(n):
    if n == 2:
        return 4 / 27, np.array([1.0, -1.0])
    points = {(k, i / 40) for k in range(n) for i in range(41)}
    free = list(range(1, n - 1))
    for _ in range(200):
        a_ub, b_ub = [], []
        for k, u in sorted(points):
            const, coef = row(n, k, u)
            fixed = const + coef[0] - coef[n - 1]  # s_0 = 1, s_{n-1} = -1
            for sign in (1, -1):
                a_ub.append(list(sign * coef[free]) + [-1.0])
                b_ub.append(-sign * fixed)
        result = linprog([0.0] * len(free) + [1.0], A_ub=a_ub, b_ub=b_ub, bounds=[(None, None)] * len(free) + [(0, None)], method="highs",
                         options={"primal_feasibility_tolerance": 1e-10, "dual_feasibility_tolerance": 1e-10})
        assert result.status == 0, result.message
        s = np.array([1.0] + list(result.x[:-1]) + [-1.0])
        m = result.x[-1]
        added = 0
        for k in range(n):
            for u in extremes(n, k, s):
                if abs(value(n, k, u, s)) > m + 1e-13 and (k, u) not in points:
                    points.add((k, u))
                    added += 1
        if not added:
            return largest(n, s), s
    raise RuntimeError("no convergence for n = %d" % n)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    max_n = int(args[0]) if args else 16
    table = {}
    print("%3s  %-14s %-14s %s" % ("n", "error factor", "limit", "weights"))
    for n in range(2, max_n + 1):
        factor, s = best(n)
        table[n] = factor
        shown = ", ".join("%.4f" % x for x in (s if n <= 9 else list(s[:4]))) + ("" if n <= 9 else ", ...")
        print("%3d  %.12f %14.4f  (%s)" % (n, factor, n**3 / factor, shown))
        # Reversing the cubic turns weights s_k into -s_{n-1-k}; the best weights
        # need not be unique, but the mean of a solution and its mirror is one
        # too, and it is antisymmetric.
        mirror = (s - s[::-1]) / 2
        assert largest(n, mirror) <= factor + 1e-12, (n, largest(n, mirror), factor)
        if n % 2 == 0:
            assert abs(factor - 4 / 27) < 1e-9, (n, factor)
        else:
            assert factor > 4 / 27
    odd = [n for n in table if n % 2]
    assert all(table[a] > table[b] for a, b in zip(odd, odd[1:])), "not decreasing over odd n"
    if "--table" in sys.argv:
        print("\nBEST_ODD_FACTOR = {")
        for n in odd:
            print("    %d: %.12f," % (n, table[n]))
        print("}")


if __name__ == "__main__":
    main()
