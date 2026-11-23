#!/usr/bin/env python3
"""Dump what fontTools' cu2qu returns for a large fixed set of inputs, one
line per input, so that two versions can be compared with diff.

Usage: PYTHONPATH=<fontTools Lib> ft_equiv.py <out-file> [cubics.json]

The inputs are generated with a fixed seed and do not depend on the
fontTools under test (cubics from fonts are read from a JSON file made
separately).
"""
import json, math, os, random, sys
from fontTools.cu2qu import cu2qu
from fontTools.cu2qu.cu2qu import curve_to_quadratic, curves_to_quadratic

out = open(sys.argv[1], 'w')
print('# cu2qu from', cu2qu.__file__, 'compiled:', cu2qu.COMPILED, file=sys.stderr)

def one(curve, tol, all_quadratic=True):
    try:
        r = curve_to_quadratic(curve, tol, all_quadratic)
        return repr(r)
    except Exception as ex:
        return type(ex).__name__

def many(curves, tols, all_quadratic=True):
    try:
        r = curves_to_quadratic(curves, tols, all_quadratic)
        return repr(r)
    except Exception as ex:
        return type(ex).__name__

rnd = random.Random(12345)
count = 0
def emit(tag, s):
    global count
    count += 1
    out.write(f'{tag} {s}\n')

# 1. Random float cubics over scales and tolerances.
for i in range(120000):
    scale = 10 ** rnd.uniform(0, 4)
    c = [(rnd.uniform(-scale, scale), rnd.uniform(-scale, scale)) for _ in range(4)]
    tol = rnd.choice([0.05, 0.25, 0.5, 1.0, 2.0, 2.048])
    emit('rf', one(c, tol))

# 2. Random integer cubics (what fonts have), smaller shapes.
for i in range(120000):
    size = rnd.choice([20, 60, 200, 600, 2048])
    x0, y0 = rnd.randint(-1000, 1000), rnd.randint(-1000, 1000)
    c = [(x0 + rnd.randint(-size, size), y0 + rnd.randint(-size, size)) for _ in range(4)]
    tol = rnd.choice([0.5, 1.0, 2.048])
    emit('ri', one(c, tol))

# 3. Third difference exactly at, or a hair off, a limit of the closed
#    form, for n = 2 (54 tol) and n = 3 (128 tol); and near the others.
def Phi(n):
    if n == 2: return 4 / 27
    g = (n - 2) / (4 * (n - 1))
    first = (4 / 27) * (1 + g) ** 3
    w2 = (0.25 - 0.5 / (n - 1)) / 3
    return max(first, (n - 3) / (4 * (n - 1)) + 2 * w2 * math.sqrt(w2))
# TIE_INPUTS=<file> also writes the inputs of this class, one JSON line each.
tie_inputs = open(os.environ['TIE_INPUTS'], 'w') if os.environ.get('TIE_INPUTS') else None
for i in range(40000):
    tol = rnd.choice([0.5, 1.0, 2.0])
    n = rnd.choice([2, 2, 3, 3, 4, 5, 6, 8, 12, 20, 40, 99, 100])
    limit = tol * n ** 3 / Phi(n)
    mag = limit * rnd.choice([1.0, 1.0, 1 - 1e-15, 1 + 1e-15, 1 - 1e-12, 1 + 1e-12, 1 - 1e-9, 1 + 1e-9,
                              1 - 1e-6, 1 + 1e-6, 0.999, 1.001])
    d = rnd.choice([(mag, 0.0), (0.0, mag), (-mag, 0.0), (0.0, -mag),
                    (mag * 0.6, mag * 0.8), (-mag * 0.8, mag * 0.6)])
    span = max(50, int(limit))
    c0 = (rnd.randint(-500, 500), rnd.randint(-500, 500))
    c1 = (c0[0] + rnd.randint(-span, span), c0[1] + rnd.randint(-span, span))
    c2 = (c1[0] + rnd.randint(-span, span), c1[1] + rnd.randint(-span, span))
    c3 = (c0[0] - 3 * c1[0] + 3 * c2[0] + d[0], c0[1] - 3 * c1[1] + 3 * c2[1] + d[1])
    emit('tie', one([c0, c1, c2, c3], tol))
    if tie_inputs: tie_inputs.write(json.dumps([[c0, c1, c2, c3], tol, n]) + '\n')

# 4. Odd inputs.
big = 1e15
odd = [
    [(0, 0), (0, 0), (0, 0), (0, 0)],
    [(5, 5), (5, 5), (5, 5), (5, 5.1)],
    [(0, 0), (1, 0), (2, 0), (3, 0)],
    [(0, 0), (3, 0), (-1, 0), (2, 0)],
    [(0, 0), (float('nan'), 0), (1, 1), (2, 2)],
    [(0, 0), (float('inf'), 0), (1, 1), (2, 2)],
    [(0, 0), (1e300, 0), (0, 1e300), (2, 2)],
    [(big, big), (big + 100, big), (big + 200, big + 300), (big + 400, big + 400)],
    [(1e9, 1e9), (1e9 + 100, 1e9), (1e9 + 200, 1e9 + 300), (1e9 + 400, 1e9 + 400)],
    [(0, 0), (1e6, 3e5), (1e6, -7e5), (3, 3e6)],
]
# (The compiled module crashes on some of these: SKIP_ODD=1 leaves them out,
# and ft_odd.py runs them one per process.)
if not os.environ.get('SKIP_ODD'):
  for c in odd:
    for tol in (1.0, 0.0, -1.0, 1e-9, 1e9, float('nan'), float('inf')):
        emit('odd', one(c, tol))
        emit('oddF', one(c, tol, False))

# 5. Several curves at once, as for compatible masters; and all_quadratic=False.
for i in range(40000):
    k = rnd.randint(1, 4)
    size = rnd.choice([30, 100, 400, 2048])
    curves = [[(rnd.randint(-size, size), rnd.randint(-size, size)) for _ in range(4)] for _ in range(k)]
    tols = [rnd.choice([0.5, 1.0, 2.048, 0.05])] * k
    if i % 5 == 0:
        tols = [rnd.choice([0.05, 0.5, 1.0, 5.0]) for _ in range(k)]
    emit('multi', many(curves, tols))
    if i % 4 == 0:
        emit('multiF', many(curves, tols, False))
        emit('oneF', one(curves[0], tols[0], False))

# 6. Cubics from fonts.
if len(sys.argv) > 2:
    cubics = json.load(open(sys.argv[2]))
    for c in cubics:
        for tol in (1.0, 0.5):
            emit('font', one(c, tol))
    # compatible pairs: consecutive cubics as if they were two masters
    for a, b in zip(cubics[0:40000:2], cubics[1:40000:2]):
        emit('fontmulti', many([a, b], [1.0, 1.0]))

out.close()
print('# inputs:', count, file=sys.stderr)
