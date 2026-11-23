#!/usr/bin/env python3
"""Scan simple cubics for ones on which fontTools' fit test for the single
quadratic recurses without end, with the fontTools on PYTHONPATH (pure
Python: a RecursionError is the sign).

Usage: PYTHONPATH=<fontTools Lib> ft_endless_scan.py FAMILY [max]
  three   (0,0) (0,0) (0,0) (X,0), X = 1..max: the tolerance is 4/27 X, and
          the doubles around it
  handle  (-3p,-3) (-2p,-2) (0,5) (0,9), p = 2..max: tangents meet at (0,0);
          the error curve has e1 = (p,1), e2 = (0,-2) and is farthest, 4p/9,
          at t = 1/3, on the x axis
"""
import math, sys

sys.setrecursionlimit(220)
from fontTools.cu2qu import cu2qu
from fontTools.cu2qu.cu2qu import curve_to_quadratic

assert not cu2qu.COMPILED, "run this on pure Python"
family = sys.argv[1]
top = int(sys.argv[2]) if len(sys.argv) > 2 else 3000


def nudge(x, k):
    for _ in range(abs(k)):
        x = math.nextafter(x, math.inf if k > 0 else -math.inf)
    return x


def outcome(cubic, tol):
    try:
        return len(curve_to_quadratic(cubic, tol)) - 2
    except RecursionError:
        return "endless"


found = 0
tests = 0
for v in range(1, top + 1):
    if family == "three":
        cubic = [(0, 0), (0, 0), (0, 0), (v, 0)]
        centre = 4 * v / 27
    else:
        if v < 2:
            continue
        cubic = [(-3 * v, -3), (-2 * v, -2), (0, 5), (0, 9)]
        centre = 4 * v / 9
    for k in range(-3, 4):
        tol = nudge(centre, k)
        tests += 1
        r = outcome(cubic, tol)
        if r == "endless":
            found += 1
            if found <= 25:
                print(f"endless: {cubic} tolerance {tol!r} (= {'4*%d/27' % v if family == 'three' else '4*%d/9' % v} {k:+d} ulp; {tol.hex()})")
print(f"{family}: {tests} tests, {found} endless")
