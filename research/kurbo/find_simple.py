#!/usr/bin/env python3
"""Find simple non-integer cubics on which upstream kurbo and fontTools
give different bits (for a regression test).

Usage: PYTHONPATH=<fontTools Lib> find_simple.py <kdiff-old binary>
"""
import random, struct, subprocess, sys
from fontTools.cu2qu.cu2qu import curve_to_quadratic


def bits(x):
    return struct.pack(">d", float(x)).hex()


def unbits(s):
    return struct.unpack(">d", bytes.fromhex(s))[0]


rnd = random.Random(5)
cases = []
for i in range(30000):
    den = rnd.choice([2, 4, 5, 10])
    size = rnd.choice([20, 50, 100])
    c = [(rnd.randint(-size * den, size * den) / den, rnd.randint(-size * den, size * den) / den) for _ in range(4)]
    cases.append((c, 1.0))
text = "".join("S %s %s\n" % (bits(tol), " ".join(bits(v) for p in c for v in p)) for c, tol in cases)
r = subprocess.run([sys.argv[1]], input=text, capture_output=True, text=True)
assert r.returncode == 0, r.stderr[-200:]
found = []
for (c, tol), line in zip(cases, r.stdout.splitlines()):
    py = curve_to_quadratic(c, tol)
    ku = [unbits(s) for s in line.split()]
    ku = list(zip(ku[0::2], ku[1::2]))
    if ku != py and len(ku) == len(py):
        digits = sum(len(repr(v)) for p in c for v in p)
        found.append((len(py), digits, c, py, ku))
found.sort(key=lambda t: (t[0], t[1]))
print("differing:", len(found), "of", len(cases))
for n, digits, c, py, ku in found[:3]:
    print("cubic", c, "->", n - 2, "quadratics")
    for a, b in zip(py, ku):
        print("   fontTools", a, "  kurbo", b, "" if a == b else "  <-- differs")
