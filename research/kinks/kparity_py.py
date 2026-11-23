#!/usr/bin/env python3
"""Do kurbo and fontTools convert real masters together to the same points?

  kparity_py.py gen <inputs> <limit> DATA.npz ...  the cubics of the masters
                                                    of stage 1, as groups
  kparity_py.py run <inputs> <outputs>              convert with the fontTools
                                                    on PYTHONPATH

Formats, with every float as the 16 hex digits of its bits (those of the
earlier kurbo-vs-fontTools harness):
  inputs   M <k> <tol> <8k coordinates>       k cubics converted together
  outputs  - if no approximation; else the points' coordinates, splines of a
           group separated by " | "
"""
import json
import random
import struct
import sys


def bits(x):
    return struct.pack(">d", x).hex()


def unbits(s):
    return struct.unpack(">d", bytes.fromhex(s))[0]


def gen(path, limit, *datasets):
    import numpy as np

    limit = int(limit)
    rnd = random.Random(7)
    out = open(path, "w")
    count = 0
    for ds in datasets:
        data = np.load(ds)
        with open(ds + ".meta.json") as f:
            tol = json.load(f)["tol"]
        pts = data["pts"]
        idx = [int(i) for i in np.nonzero(data["kind"] == 1)[0]]
        if len(idx) > limit:
            idx = sorted(rnd.sample(idx, limit))
        M = pts.shape[1]
        for s in idx:
            coords = []
            for m in range(M):
                for k in range(4):
                    z = pts[s, m, k]
                    coords.append(bits(float(z.real)))
                    coords.append(bits(float(z.imag)))
            out.write("M %d %s %s\n" % (M, bits(float(tol)), " ".join(coords)))
            count += 1
    out.close()
    print("groups:", count, file=sys.stderr)


def run(inputs, outputs):
    from fontTools.cu2qu import cu2qu
    from fontTools.cu2qu.cu2qu import curves_to_quadratic

    print("cu2qu from", cu2qu.__file__, "compiled:", cu2qu.COMPILED, file=sys.stderr)
    out = open(outputs, "w")

    def fmt(spline):
        return " ".join(bits(v) for p in spline for v in p)

    for line in open(inputs):
        f = line.split()
        k = int(f[1])
        v = [unbits(s) for s in f[2:]]
        curves = [
            [(v[1 + 8 * i + 2 * j], v[2 + 8 * i + 2 * j]) for j in range(4)]
            for i in range(k)
        ]
        try:
            out.write(
                " | ".join(fmt(s) for s in curves_to_quadratic(curves, [v[0]] * k))
                + "\n"
            )
        except Exception:
            out.write("-\n")
    out.close()


if __name__ == "__main__":
    if sys.argv[1] == "gen":
        gen(*sys.argv[2:])
    else:
        run(*sys.argv[2:])
