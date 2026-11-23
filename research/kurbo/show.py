#!/usr/bin/env python3
"""Show one input of the kurbo/fontTools comparison and where two outputs
differ, and tabulate over all differing lines of a class by number of
quadratics and by which point indices differ.

Usage: show.py inputs.txt a.txt b.txt <line number> [first last]
"""
import collections, struct, sys


def unbits(s):
    return struct.unpack(">d", bytes.fromhex(s))[0]


def parse(line):
    line = line.strip()
    if line == "-":
        return None
    return [unbits(s) for s in line.split(" | ")[0].split()]


inputs = open(sys.argv[1]).read().splitlines()
a = open(sys.argv[2]).read().splitlines()
b = open(sys.argv[3]).read().splitlines()
i = int(sys.argv[4]) - 1
f = inputs[i].split()
v = [unbits(s) for s in f[1:]]
print("accuracy", repr(v[0]), "cubic", [(v[1], v[2]), (v[3], v[4]), (v[5], v[6]), (v[7], v[8])])
pa, pb = parse(a[i]), parse(b[i])
print("quadratics:", len(pa) // 2 - 2, len(pb) // 2 - 2)
for j in range(0, len(pa), 2):
    mark = "" if (pa[j], pa[j + 1]) == (pb[j], pb[j + 1]) else "   <-- differs"
    print(f"  point {j // 2}: {pa[j]!r}, {pa[j + 1]!r}   |   {pb[j]!r}, {pb[j + 1]!r}{mark}")

if len(sys.argv) > 6:
    first, last = int(sys.argv[5]) - 1, int(sys.argv[6])
    by_n = collections.Counter()
    all_n = collections.Counter()
    for k in range(first, min(last, len(a), len(b))):
        pa, pb = parse(a[k]), parse(b[k])
        if pa is None or pb is None or len(pa) != len(pb):
            continue
        n = len(pa) // 2 - 2
        all_n[n] += 1
        if pa != pb:
            by_n[n] += 1
    print("lines by number of quadratics: differing / all")
    for n in sorted(all_n):
        print(f"  n={n}: {by_n[n]} / {all_n[n]}")
