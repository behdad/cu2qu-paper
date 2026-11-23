#!/usr/bin/env python3
"""Print one cubic and what a method made of it.

  dump_one.py <method.tol> <index | worst | worst:k>
"""
import json
import os
import struct
import sys

S = "/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad"
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


name, which = sys.argv[1], sys.argv[2]
cubics = json.load(open(os.path.join(S, "ft_cubics.json")))
res = read_results(os.path.join(B, "res", name + ".bin"))
per = [l.split() for l in open(os.path.join(B, "percubic", name + ".txt"))][1:]
if which.startswith("worst"):
    k = int(which.split(":")[1]) if ":" in which else 0
    order = sorted(range(len(per)), key=lambda i: -max(float(per[i][4]), float(per[i][5])))
    idxs = [order[k]]
else:
    idxs = [int(which)]
for i in idxs:
    print("cubic", i, cubics[i])
    print("  nseg %s haus_c2s %s haus_s2c %s eq %s tan %s %s kink %s bbox %s d3 %s" % tuple(per[i][j] for j in (1, 4, 5, 6, 7, 8, 9, 10, 11)))
    for s in res[i]:
        kind = "line" if s[0] == 1 else "quad"
        print("  %s (%.4f, %.4f) ctrl (%.4f, %.4f) -> (%.4f, %.4f)" % ((kind,) + tuple(s[1:])))
