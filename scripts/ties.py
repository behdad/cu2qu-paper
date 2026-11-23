#!/usr/bin/env python3
"""What cu2qu's recursive test does with a spline exactly at its limit.

Usage: ties.py REPO REV_BEFORE REV_AFTER

Random cubics with integer coordinates whose third difference is exactly
L(n) times an exactly representable tolerance, for the three limits of
Table 1 that are binary fractions: 54 (two quadratics), 128 (three) and
843.75 (six). The spline of n quadratics touches the tolerance without
crossing it, so n is the answer. The revision before (pure Python, read out
of the git repository REPO) searches with the recursive test, which works in
floating point on halves of halves, and does one of three things: accepts
the spline, rejects it and goes on to n + 1, or does not return (a
RecursionError here; a crash when compiled). The revision after computes n.
"""

import math
import random
import sys

import compare_revs

# directions with integer components for an integer or half-integer length
TRIPLES = [(1, 0, 1), (0, 1, 1), (3, 4, 5), (4, 3, 5), (-3, 4, 5), (3, -4, 5), (5, 12, 13), (12, 5, 13), (8, 15, 17), (7, 24, 25)]
LIMITS = [(2, 54.0), (3, 128.0), (6, 843.75)]
TOLERANCES = [0.5, 1.0, 2.0, 4.0, 8.0]
PER_DIRECTION = 300


def outcome(module, cubic, tolerance):
    try:
        return len(module.curve_to_quadratic(cubic, tolerance)) - 2
    except RecursionError:
        return None


def main():
    repo, rev_a, rev_b = sys.argv[1:4]
    before, _ = compare_revs.load_revision(repo, rev_a, "cu2qu_before")
    after, _ = compare_revs.load_revision(repo, rev_b, "cu2qu_after")
    rnd = random.Random(1)
    print("before: %s (pure Python)\nafter:  %s (pure Python)\n" % (rev_a, rev_b))
    print("%-22s %9s %8s | %9s %9s %14s | %s" % ("third difference", "tolerance", "cubics", "accepted", "rejected", "did not return", "after"))
    for n, limit in LIMITS:
        totals = [0, 0, 0, 0]
        for tolerance in TOLERANCES:
            length = limit * tolerance
            tally = [0, 0, 0, 0]
            agree = True
            for a, b, h in TRIPLES:
                dx, dy = length * a / h, length * b / h
                if dx != int(dx) or dy != int(dy):
                    continue
                dx, dy = int(dx), int(dy)
                assert math.hypot(dx, dy) == length
                for _ in range(PER_DIRECTION):
                    c0, c1, c2 = [(rnd.randint(-500, 500), rnd.randint(-500, 500)) for _ in range(3)]
                    c3 = (dx + 3 * c2[0] - 3 * c1[0] + c0[0], dy + 3 * c2[1] - 3 * c1[1] + c0[1])
                    cubic = [c0, c1, c2, c3]
                    tally[0] += 1
                    result = outcome(before, cubic, tolerance)
                    if result is None:
                        tally[3] += 1
                    elif result == n:
                        tally[1] += 1
                    else:
                        assert result == n + 1, (cubic, tolerance, result)
                        tally[2] += 1
                    agree = agree and outcome(after, cubic, tolerance) == n
            for i in range(4):
                totals[i] += tally[i]
            if tally[0]:
                print("%-22s %9g %8d | %9d %9d %14d | %s" % ("%g x tolerance" % limit, tolerance, tally[0], tally[1], tally[2], tally[3], "all %d quadratics" % n if agree else "DIFFERS"))
        print("%-22s %9s %8d | %9d %9d %14d |" % ("%g x tolerance" % limit, "all", *totals))
        print()


if __name__ == "__main__":
    main()
