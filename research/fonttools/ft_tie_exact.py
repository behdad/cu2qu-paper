#!/usr/bin/env python3
"""For the inputs where old and new cu2qu disagree, how far is the cubic's
third difference from the limit of the closed form, computed exactly?

Usage: ft_tie_exact.py old.txt new.txt tie_inputs.jsonl cubics.json

Every differing result should be a cubic whose |d3| is within rounding of
tol * n^3 / Phi(n) for the n in question; this prints the largest relative
distance from the limit among them, by kind of difference, and whether the
new code's answer is the exactly right one.
"""
import collections, json, math, sys
from decimal import Decimal, getcontext
from fractions import Fraction

getcontext().prec = 60


def Phi(n):
    g = (n - 2) / (4 * (n - 1))
    error = (4 / 27) * (1 + g) * (1 + g) * (1 + g)
    if n > 3:
        c = (n - 3) / (4 * (n - 1))
        w = c / 3
        error = max(error, c + 2 * w * math.sqrt(w))
    return error


def Phi_exact(n):
    """The limit's factor in 60 digits."""
    n = Decimal(n)
    g = (n - 2) / (4 * (n - 1))
    first = Decimal(4) / 27 * (1 + g) ** 3
    if n > 3:
        c = (n - 3) / (4 * (n - 1))
        w = c / 3
        return max(first, c + 2 * w * w.sqrt())
    return first


def rel_excess(cubic, tol, n):
    """|d3| / (tol n^3 / Phi(n)) - 1, exactly (60 digits)."""
    (x0, y0), (x1, y1), (x2, y2), (x3, y3) = [(Fraction(x), Fraction(y)) for x, y in cubic]
    dx = x3 - 3 * x2 + 3 * x1 - x0
    dy = y3 - 3 * y2 + 3 * y1 - y0
    d2 = dx * dx + dy * dy
    d3 = (Decimal(d2.numerator) / Decimal(d2.denominator)).sqrt()
    limit = Decimal(Fraction(tol).numerator) / Decimal(Fraction(tol).denominator) * Decimal(n) ** 3 / Phi_exact(n)
    return d3 / limit - 1


def count(line):
    tag, r = line.split(' ', 1)
    if r.startswith('['):
        return tag, r.count('(')  # points
    return tag, r.strip()


old = open(sys.argv[1]).read().splitlines()
new = open(sys.argv[2]).read().splitlines()
ties = [json.loads(l) for l in open(sys.argv[3])]
cubics = json.load(open(sys.argv[4]))
assert len(old) == len(new)

# line number -> input, for the classes that can differ
inputs = {}
tie_i = font_i = 0
for i, line in enumerate(old):
    tag = line.split(' ', 1)[0]
    if tag == 'tie':
        inputs[i] = (ties[tie_i][0], ties[tie_i][1])
        tie_i += 1
    elif tag == 'font':
        inputs[i] = (cubics[font_i // 2], (1.0, 0.5)[font_i % 2])
        font_i += 1
assert tie_i == len(ties) and font_i == 2 * len(cubics), (tie_i, font_i)

MAX_N = 100
stats = collections.defaultdict(lambda: [0, Decimal(0), 0, 0, Decimal(0)])  # count, max |rel|, new right, exact ties, max |rel| where new is not right
other = collections.Counter()
for i, (a, b) in enumerate(zip(old, new)):
    if a == b:
        continue
    tag, ra = count(a)
    _, rb = count(b)
    if i not in inputs:
        other[tag] += 1
        continue
    cubic, tol = inputs[i]
    na = ra - 2 if isinstance(ra, int) else MAX_N + 1   # ApproxNotFoundError: more than MAX_N
    nb = rb - 2 if isinstance(rb, int) else MAX_N + 1
    if not (isinstance(ra, int) or ra == 'ApproxNotFoundError') or not (isinstance(rb, int) or rb == 'ApproxNotFoundError'):
        other[f'{tag}: {ra} -> {rb}'] += 1
        continue
    if abs(na - nb) != 1:
        other[f'{tag}: n {na} -> {nb}'] += 1
        continue
    n = min(na, nb)                     # the n in dispute: does the spline of n fit?
    rel = rel_excess(cubic, tol, n)
    kind = (tag, 'new takes one fewer' if nb < na else 'new takes one more')
    st = stats[kind]
    st[0] += 1
    st[1] = max(st[1], abs(rel))
    fits = rel <= 0
    st[2] += (fits == (nb == n))        # is the new code's answer the exact one?
    st[3] += (rel == 0)
    if fits != (nb == n):
        st[4] = max(st[4], abs(rel))
print('differing inputs outside the two classes, or of another kind:', dict(other))
print(f"{'class':6s} {'difference':22s} {'count':>6s} {'largest |d3/limit - 1|':>24s} {'new exactly right':>18s} {'exact ties':>11s} {'largest where new is not':>25s}")
for (tag, kind), (c, m, right, ties_, mw) in sorted(stats.items()):
    print(f'{tag:6s} {kind:22s} {c:6d} {float(m):24.3e} {right:18d} {ties_:11d} {float(mw):25.3e}')
