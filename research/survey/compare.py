#!/usr/bin/env python3
"""Pairwise comparison of results files: are two implementations' outputs the same?

  compare.py <a> <b> [tol ...]     e.g. compare.py fonttools_cu2qu kurbo_approx_spline
"""
import json
import os
import struct
import sys

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


def main():
    a, b = sys.argv[1], sys.argv[2]
    tols = sys.argv[3:] or ["1.0", "0.5", "0.1"]
    cubics = json.load(open(os.path.join(S, "ft_cubics.json")))
    for tol in tols:
        ra = read_results(os.path.join(B, "res", f"{a}.{tol}.bin"))
        rb = read_results(os.path.join(B, "res", f"{b}.{tol}.bin"))
        same_bits = 0
        same_n = 0
        maxdev = 0.0
        ndiff = []
        devs = []
        for i, (x, y) in enumerate(zip(ra, rb)):
            if x == y:
                same_bits += 1
                same_n += 1
                continue
            if len(x) != len(y):
                ndiff.append((i, len(x), len(y)))
                continue
            same_n += 1
            d = max(abs(p - q) for s, t in zip(x, y) for p, q in zip(s[1:], t[1:]))
            maxdev = max(maxdev, d)
            devs.append((d, i))
        devs.sort(reverse=True)
        print(
            f"tol {tol}: {a} vs {b}: bit-identical {same_bits}/{len(ra)}, same count {same_n}, "
            f"different count {len(ndiff)}, largest coordinate difference among same-count {maxdev:.3e}"
        )
        for i, na, nb in ndiff[:12]:
            c = cubics[i]
            d3 = (
                (c[3][0] - 3 * c[2][0] + 3 * c[1][0] - c[0][0]) ** 2
                + (c[3][1] - 3 * c[2][1] + 3 * c[1][1] - c[0][1]) ** 2
            ) ** 0.5
            print(f"    cubic {i} {c}: {na} vs {nb} segments; |d3|/tol = {d3 / float(tol):.6f}")
        if len(ndiff) > 12:
            print(f"    ... and {len(ndiff) - 12} more")
        for d, i in devs[:5]:
            print(f"    cubic {i} {cubics[i]}: same count {len(ra[i])}, coordinates differ by {d:.3e}")


main()
