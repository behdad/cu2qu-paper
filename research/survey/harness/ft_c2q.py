#!/usr/bin/env python3
"""Run upstream fontTools cu2qu (pure Python, from the clean export) over the corpus.

  PYTHONPATH=<S>/ft-old/Lib python3 ft_c2q.py <corpus.json> <tol> <out.bin>

Output: "C2QRES01", u32 count, per cubic u32 nseg, nseg * {u32 kind=2, 6 f64}.
The spline returned by curve_to_quadratic (first point, off-curve points, last
point) is expanded to explicit quadratics with the implied on-curve points.
"""
import json
import struct
import sys
import time

import fontTools.cu2qu.cu2qu as cu2qu
from fontTools.cu2qu.cu2qu import curve_to_quadratic


def main():
    corpus, tol, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    cubics = json.load(open(corpus))
    assert not cu2qu.COMPILED, "expected the pure-Python module"
    t0 = time.perf_counter()
    splines = [curve_to_quadratic(c, tol) for c in cubics]
    dt = time.perf_counter() - t0
    buf = bytearray()
    buf += b"C2QRES01"
    buf += struct.pack("<I", len(cubics))
    nq = 0
    for sp in splines:
        n = len(sp) - 2
        nq += n
        buf += struct.pack("<I", n)
        prev = sp[0]
        for k in range(n):
            off = sp[k + 1]
            if k + 1 < n:
                nxt = sp[k + 2]
                on = ((off[0] + nxt[0]) * 0.5, (off[1] + nxt[1]) * 0.5)
            else:
                on = sp[-1]
            buf += struct.pack("<I6d", 2, prev[0], prev[1], off[0], off[1], on[0], on[1])
            prev = on
    open(out, "wb").write(buf)
    print(
        json.dumps(
            {
                "method": "fonttools_cu2qu_purepython",
                "tol": tol,
                "ns_per_cubic": round(dt * 1e9 / len(cubics), 1),
                "quads": nq,
                "path": cu2qu.__file__,
            }
        )
    )


main()
