#!/usr/bin/env python3
"""Independent check of the evaluator's Hausdorff distance: brute force with numpy.

  verify_haus.py <method.tol> [count] [seed]

Samples each curve at 3001 points per segment and takes point-to-polyline
distances both ways (no root finding, no refinement), then compares with the
harness's per-cubic output.
"""
import json
import os
import random
import struct
import sys

import numpy as np

S = "/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad"
B = "/home/behdad/hb/cu2qu-work/research-survey"


def read_results(path):
    buf = open(path, "rb").read()
    assert buf[:8] == b"C2QRES01"
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


def cubic_pts(c, m):
    t = np.linspace(0, 1, m)[:, None]
    p = np.array(c, dtype=float)
    return (1 - t) ** 3 * p[0] + 3 * (1 - t) ** 2 * t * p[1] + 3 * (1 - t) * t**2 * p[2] + t**3 * p[3]


def seg_pts(s, m):
    t = np.linspace(0, 1, m)[:, None]
    p0 = np.array(s[1:3])
    p1 = np.array(s[3:5])
    p2 = np.array(s[5:7])
    if s[0] == 1:
        return (1 - t) * p0 + t * p2
    return (1 - t) ** 2 * p0 + 2 * (1 - t) * t * p1 + t**2 * p2


def pt_polyline_dist(P, Q):
    """max over points P of the distance to polyline Q."""
    a = Q[:-1]
    b = Q[1:]
    ab = b - a
    ab2 = (ab * ab).sum(1)
    ab2[ab2 == 0] = 1.0
    worst = 0.0
    for i in range(0, len(P), 256):
        p = P[i : i + 256, None, :]
        u = (((p - a) * ab).sum(2) / ab2).clip(0, 1)
        d = p - (a + u[:, :, None] * ab)
        worst = max(worst, np.sqrt((d * d).sum(2).min(1)).max())
    return worst


def main():
    name = sys.argv[1]
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    cubics = json.load(open(os.path.join(S, "ft_cubics.json")))
    res = read_results(os.path.join(B, "res", name + ".bin"))
    per = [l.split() for l in open(os.path.join(B, "percubic", name + ".txt"))][1:]
    rnd = random.Random(seed)
    idxs = rnd.sample(range(len(cubics)), count)
    worst_rel = 0.0
    worst_abs = 0.0
    for i in idxs:
        segs = res[i]
        C = cubic_pts(cubics[i], 3001 * max(1, len(segs)) // 2)
        Sp = np.concatenate([seg_pts(s, 3001) for s in segs])
        # polyline through the spline's samples: drop the jump-free duplicates only
        h1 = pt_polyline_dist(C, Sp)
        h2 = pt_polyline_dist(Sp, C)
        mine_c2s, mine_s2c = float(per[i][4]), float(per[i][5])
        for a, b in ((h1, mine_c2s), (h2, mine_s2c)):
            worst_abs = max(worst_abs, abs(a - b))
            if max(a, b) > 1e-6:
                worst_rel = max(worst_rel, abs(a - b) / max(a, b))
    print(f"{name}: {count} cubics, largest |brute - harness| = {worst_abs:.3e} units, largest relative = {worst_rel:.3e}")


main()
