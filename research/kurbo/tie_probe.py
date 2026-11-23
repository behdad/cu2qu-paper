#!/usr/bin/env python3
"""Does a kdiff binary survive cubics whose third difference is exactly at a
limit of the closed form? For each (third difference, accuracy) it makes
random integer cubics with exactly that third difference and runs each in
its own process.

Usage: tie_probe.py <binary> [count per case]
"""
import random, struct, subprocess, sys

binary = sys.argv[1]
count = int(sys.argv[2]) if len(sys.argv) > 2 else 200


def bits(x):
    return struct.pack(">d", float(x)).hex()


cases = [
    # (third difference, accuracy, what)
    ((54, 0), 1.0, "n=2 limit, 54 x 1"),
    ((128, 0), 1.0, "n=3 limit, 128 x 1"),
    ((0, 128), 1.0, "n=3 limit, vertical"),
    ((27, 0), 0.5, "n=2 limit, 54 x 0.5"),
    ((64, 0), 0.5, "n=3 limit, 128 x 0.5"),
    ((162, 216), 5.0, "n=2 limit, 54 x 5, diagonal"),
    ((384, 512), 5.0, "n=3 limit, 128 x 5, diagonal"),
    ((1728, 0), 2.048, "n=6 limit, 843.75 x 2.048 (2048 upm default)"),
    ((0, 1728), 2.048, "n=6 limit, vertical"),
    ((675, 0), 0.8, "n=6 limit, 843.75 x 0.8"),
    ((405, 540), 0.8, "n=6 limit, 843.75 x 0.8, diagonal"),
    ((3375, 0), 4.0, "n=6 limit, 843.75 x 4"),
    ((2025, 2700), 4.0, "n=6 limit, 843.75 x 4, diagonal"),
    ((253.125, 337.5), 0.5, "n=6 limit, 843.75 x 0.5, diagonal, eighths"),
    ((421.875, 0), 0.5, "n=6 limit, 843.75 x 0.5, eighths"),
    ((843.75, 0), 1.0, "n=6 limit, 843.75 x 1, quarters"),
    ((506.25, 675), 1.0, "n=6 limit, 843.75 x 1, diagonal, quarters"),
]
rnd = random.Random(7)
for d3, acc, what in cases:
    crashes = 0
    first = None
    results = {}
    for i in range(count):
        size = rnd.choice([200, 600, 2000])
        c0 = (rnd.randint(-500, 500), rnd.randint(-500, 500))
        c1 = (c0[0] + rnd.randint(-size, size), c0[1] + rnd.randint(-size, size))
        c2 = (c1[0] + rnd.randint(-size, size), c1[1] + rnd.randint(-size, size))
        c3 = (c0[0] - 3 * c1[0] + 3 * c2[0] + d3[0], c0[1] - 3 * c1[1] + 3 * c2[1] + d3[1])
        line = "S %s %s\n" % (bits(acc), " ".join(bits(v) for p in (c0, c1, c2, c3) for v in p))
        r = subprocess.run([binary], input=line, capture_output=True, text=True)
        if r.returncode != 0:
            crashes += 1
            if first is None:
                first = (c0, c1, c2, c3)
        else:
            out = r.stdout.strip()
            n = "none" if out == "-" else len(out.split()) // 2 - 2
            results[n] = results.get(n, 0) + 1
    print(f"{what:48s} crashes {crashes:3d}/{count}  quadratics {dict(sorted(results.items(), key=str))}" + (f"  first crash {first}" if first else ""))
