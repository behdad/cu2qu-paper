#!/usr/bin/env python3
"""Run single cubics through a kdiff binary, one process each.

Usage: one.py <binary> "<accuracy> x0 y0 x1 y1 x2 y2 x3 y3" ...
Prints the number of quadratics, "none", or CRASH for each.
"""
import struct, subprocess, sys


def bits(x):
    return struct.pack(">d", float(x)).hex()


for arg in sys.argv[2:]:
    v = [float(s) for s in arg.split()]
    line = "S " + " ".join(bits(x) for x in v) + "\n"
    r = subprocess.run([sys.argv[1]], input=line, capture_output=True, text=True)
    if r.returncode != 0:
        res = "CRASH (%s)" % (r.stderr.strip().splitlines()[-1] if r.stderr.strip() else r.returncode)
    else:
        out = r.stdout.strip()
        res = "none" if out == "-" else "%d quadratics" % (len(out.split()) // 2 - 2)
    print(f"{arg:60s} {res}")
