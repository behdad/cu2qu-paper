#!/usr/bin/env python3
"""Write the reproducers of the endless fit test in the comparison
harness's format ("S <tol> <8 coords>", floats as hex bits).

Usage: mk_inputs.py <out>
"""
import struct, sys


def bits(x):
    return struct.pack(">d", float(x)).hex()


cases = [
    ([(0, 0), (0, 0), (0, 0), (25, 0)], 100 / 27),
    ([(-45, -3), (-30, -2), (0, 5), (0, 9)], 20 / 3),
    ([(0, 0), (0, 0), (0, 0), (1282, 0)], 189.92592592592595),
    ([(47.09375, 0), (0, 0), (0, 0), (0, 0)], 6.9768518518518503),
    ([(1212, 1113), (0, 0), (0, 0), (0, 0)], 243.77970221057475),
    # the two of kurbo's own test
    ([(-31, 97), (-51, 93), (-91, 76), (-139, 52)], 1.987615979999813),
    ([(-20, 7), (-20, -27.5), (-12, -32.5), (4, -32)], 3.5555555555555554),
]
with open(sys.argv[1], "w") as f:
    for c, tol in cases:
        f.write("S %s %s\n" % (bits(tol), " ".join(bits(v) for p in c for v in p)))
