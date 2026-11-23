#!/usr/bin/env python3
"""How far are the interior on-curve points of a method's output from the
midpoints of the off-curve points on either side (0 = TrueType can imply them)?

  implied_dev.py <method.tol> [...]
"""
import os
import struct
import sys

import numpy as np

B = "/home/behdad/hb/cu2qu-work/research-survey"


def read_results(path):
    buf = open(path, "rb").read()
    (count,) = struct.unpack_from("<I", buf, 8)
    pos = 12
    out = []
    for _ in range(count):
        (n,) = struct.unpack_from("<I", buf, pos)
        pos += 4
        segs = []
        for _ in range(n):
            segs.append(struct.unpack_from("<I6d", buf, pos))
            pos += 52
        out.append(segs)
    return out


for name in sys.argv[1:]:
    res = read_results(os.path.join(B, "res", name + ".bin"))
    devs = []
    for segs in res:
        for a, b in zip(segs, segs[1:]):
            if a[0] != 2 or b[0] != 2:
                continue
            mx = (a[3] + b[3]) / 2
            my = (a[4] + b[4]) / 2
            devs.append(((a[5] - mx) ** 2 + (a[6] - my) ** 2) ** 0.5)
    d = np.array(devs)
    print(
        f"{name}: {len(d)} interior joints between quadratics; on-curve point is the midpoint (to 1e-6) at "
        f"{int((d <= 1e-6).sum())} ({100 * (d <= 1e-6).mean():.1f}%); distance from the midpoint: "
        f"median {np.median(d):.3f}, p90 {np.percentile(d, 90):.3f}, max {d.max():.3f} units"
    )
