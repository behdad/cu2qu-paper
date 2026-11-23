#!/usr/bin/env python3
"""What fontTools' cu2qu does on odd inputs (NaN, infinity, zero or negative
tolerance, huge coordinates), one process per input so that a crash of the
interpreter is an outcome and not the end of the run.

Usage: PYTHONPATH=<fontTools Lib> ft_odd.py [python-to-use]
Prints one line per input: the outcome, summarized.
"""
import collections, os, subprocess, sys

python = sys.argv[1] if len(sys.argv) > 1 else sys.executable
big = 1e15
odd = [
    ("all equal", "[(0, 0), (0, 0), (0, 0), (0, 0)]"),
    ("three equal", "[(5, 5), (5, 5), (5, 5), (5, 5.1)]"),
    ("straight", "[(0, 0), (1, 0), (2, 0), (3, 0)]"),
    ("straight, back and forth", "[(0, 0), (3, 0), (-1, 0), (2, 0)]"),
    ("NaN coordinate", "[(0, 0), (float('nan'), 0), (1, 1), (2, 2)]"),
    ("NaN last point", "[(0, 0), (1, 5), (1, 1), (2, float('nan'))]"),
    ("infinite coordinate", "[(0, 0), (float('inf'), 0), (1, 1), (2, 2)]"),
    ("1e300 coordinates", "[(0, 0), (1e300, 0), (0, 1e300), (2, 2)]"),
    ("1e15 offset", f"[({big}, {big}), ({big + 100}, {big}), ({big + 200}, {big + 300}), ({big + 400}, {big + 400})]"),
    ("large", "[(0, 0), (1e6, 3e5), (1e6, -7e5), (3, 3e6)]"),
]
tols = ["1.0", "0.0", "-1.0", "1e-9", "1e9", "float('nan')", "float('inf')"]

prog = r'''
import sys
from fontTools.cu2qu.cu2qu import curve_to_quadratic, curves_to_quadratic
c = %s
tol = %s
allq = %s
multi = %s
try:
    r = curves_to_quadratic([c, c], [tol, tol], allq)[0] if multi else curve_to_quadratic(c, tol, allq)
    bad = any(x != x or abs(x) == float('inf') for p in r for x in p)
    print('%%d points%%s' %% (len(r), ' (not finite)' if bad else ''))
except Exception as ex:
    print(type(ex).__name__)
'''
counts = collections.Counter()
for name, c in odd:
    for tol in tols:
        for allq in (True, False):
            for multi in (False, True):
                r = subprocess.run([python, "-c", prog % (c, tol, allq, multi)], capture_output=True, text=True)
                out = r.stdout.strip() if r.returncode == 0 else f"CRASH (exit {r.returncode})"
                counts[out.split(' (')[0] if out.startswith('CRASH') else out] += 1
                print(f"{name:26s} tol {tol:13s} all_quadratic={allq!s:5s} {'curves' if multi else 'curve '}  {out}")
print("# summary:", dict(counts))
