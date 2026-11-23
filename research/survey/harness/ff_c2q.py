#!/usr/bin/env python3
"""Run FontForge's own cubic -> quadratic conversion (splineorder2.c, through the
installed `fontforge` Python module) over the corpus, one open contour per cubic.

  python3 ff_c2q.py <corpus.json> <out.bin> [start end]

FontForge has no tolerance parameter: the test in comparespline() allows one
em-unit (1.5 when it rounds control points to integers, which it does when the
end points are integers, as in this corpus).

Output: "C2QRES01", u32 count, per cubic u32 nseg, nseg * {u32 kind, 6 f64};
kind 1 = line (two consecutive on-curve points), 2 = quadratic.
"""
import json
import struct
import sys
import time

import fontforge


def convert(c):
    ct = fontforge.contour()
    ct.moveTo(c[0][0], c[0][1])
    ct.cubicTo((c[1][0], c[1][1]), (c[2][0], c[2][1]), (c[3][0], c[3][1]))
    ct.is_quadratic = True
    return [(p.x, p.y, p.on_curve) for p in ct]


def to_segs(pts):
    segs = []
    assert pts[0][2], pts
    prev = (pts[0][0], pts[0][1])
    i = 1
    n = len(pts)
    while i < n:
        x, y, on = pts[i]
        if on:
            segs.append((1, prev[0], prev[1], (prev[0] + x) / 2, (prev[1] + y) / 2, x, y))
            prev = (x, y)
            i += 1
        else:
            assert i + 1 < n, pts
            nx, ny, non = pts[i + 1]
            if non:
                end = (nx, ny)
                i += 2
            else:  # two off-curve points in a row: implied on-curve point
                end = ((x + nx) / 2, (y + ny) / 2)
                i += 1
            segs.append((2, prev[0], prev[1], x, y, end[0], end[1]))
            prev = end
    return segs


def main():
    corpus, out = sys.argv[1], sys.argv[2]
    cubics = json.load(open(corpus))
    if len(sys.argv) > 4:
        cubics = cubics[int(sys.argv[3]) : int(sys.argv[4])]
    t0 = time.perf_counter()
    raw = [convert(c) for c in cubics]
    dt = time.perf_counter() - t0
    buf = bytearray(b"C2QRES01")
    buf += struct.pack("<I", len(cubics))
    nq = nl = 0
    moved = 0
    for c, pts in zip(cubics, raw):
        segs = to_segs(pts)
        if (pts[0][0], pts[0][1]) != tuple(c[0]) or (pts[-1][0], pts[-1][1]) != tuple(c[3]):
            moved += 1
        buf += struct.pack("<I", len(segs))
        for s in segs:
            buf += struct.pack("<I6d", *s)
            if s[0] == 2:
                nq += 1
            else:
                nl += 1
    open(out, "wb").write(buf)
    print(
        json.dumps(
            {
                "method": "fontforge_%s" % fontforge.version(),
                "ns_per_cubic_incl_python_objects": round(dt * 1e9 / len(cubics), 1),
                "quads": nq,
                "lines": nl,
                "endpoints_moved": moved,
            }
        )
    )


main()
