#!/usr/bin/env python3
"""Classify the differences between two output files of the kurbo/fontTools
comparison (see kdiff_py.py), by input class.

Usage: analyze.py inputs.txt a.txt b.txt [max_lines]
"""
import collections, struct, sys


def unbits(s):
    return struct.unpack(">d", bytes.fromhex(s))[0]


def parse(line):
    line = line.strip()
    if line == "-":
        return None
    return [[unbits(s) for s in part.split()] for part in line.split(" | ")]


inputs = open(sys.argv[1])
a = open(sys.argv[2])
b = open(sys.argv[3])
limit = int(sys.argv[4]) if len(sys.argv) > 4 else None

# classes by line number, as kdiff_py.py gen writes them
bounds = [(120000, "random float"), (240000, "random int"), (280000, "at a limit"), (340000, "nearly quadratic"), (380000, "groups"), (535066, "font"), (555066, "font pairs")]


def cls(i):
    for end, name in bounds:
        if i < end:
            return name
    return "?"


kinds = collections.defaultdict(collections.Counter)
worst = collections.defaultdict(float)
examples = {}
n = 0
for i, (la, lb) in enumerate(zip(a, b)):
    if limit is not None and i >= limit:
        break
    n += 1
    c = cls(i)
    if la == lb:
        kinds[c]["identical"] += 1
        continue
    pa, pb = parse(la), parse(lb)
    if pa is None or pb is None:
        kind = "one has no result"
    elif [len(s) for s in pa] != [len(s) for s in pb]:
        da = len(pa[0]) // 2 - 2
        db = len(pb[0]) // 2 - 2
        kind = "b has one quadratic fewer" if db == da - 1 else "b has one quadratic more" if db == da + 1 else f"counts {da} vs {db}"
    else:
        d = max(abs(x - y) for sa, sb in zip(pa, pb) for x, y in zip(sa, sb))
        m = max(max(abs(x) for s in pa for x in s), 1.0)
        rel = d / m
        worst[c] = max(worst[c], rel)
        kind = "same count, zero sign only" if d == 0 else "same count, within 1e-12" if rel <= 1e-12 else "same count, within 1e-6" if rel <= 1e-6 else "same count, LARGER"
    kinds[c][kind] += 1
    examples.setdefault((c, kind), i)
print("lines compared:", n)
for c in [name for _, name in bounds]:
    if c in kinds:
        k = kinds[c]
        total = sum(k.values())
        rest = ", ".join(f"{name}: {v}" for name, v in sorted(k.items()) if name != "identical")
        print(f"{c:18s} {total:7d}  identical {k['identical']:7d}  {rest}   (largest relative coordinate difference {worst[c]:.2e})")
print("first example line of each kind:", {f"{c}/{k}": i + 1 for (c, k), i in examples.items()})
