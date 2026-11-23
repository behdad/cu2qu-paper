#!/usr/bin/env python3
"""What a glyf table would hold: round the off-curve points of an implied-on-curve
spline to integers (round half up, as fontTools' otRound does) and re-derive the
implied on-curve points from the rounded ones.  Only for results whose interior
on-curve points are all implied.

  round_offcurves.py <in.bin> <out.bin>
"""
import math
import struct
import sys


def ot_round(v):
    return float(math.floor(v + 0.5))


buf = open(sys.argv[1], "rb").read()
assert buf[:8] == b"C2QRES01"
(count,) = struct.unpack_from("<I", buf, 8)
pos = 12
out = bytearray(b"C2QRES01")
out += struct.pack("<I", count)
moved = 0
for _ in range(count):
    (n,) = struct.unpack_from("<I", buf, pos)
    pos += 4
    segs = []
    for _ in range(n):
        segs.append(struct.unpack_from("<I6d", buf, pos))
        pos += 52
    offs = [(ot_round(s[3]), ot_round(s[4])) for s in segs]
    p0 = (segs[0][1], segs[0][2])
    p3 = (segs[-1][5], segs[-1][6])
    out += struct.pack("<I", n)
    prev = p0
    for k in range(n):
        if k + 1 < n:
            on = ((offs[k][0] + offs[k + 1][0]) / 2, (offs[k][1] + offs[k + 1][1]) / 2)
        else:
            on = p3
        out += struct.pack("<I6d", 2, prev[0], prev[1], offs[k][0], offs[k][1], on[0], on[1])
        prev = on
open(sys.argv[2], "wb").write(out)
