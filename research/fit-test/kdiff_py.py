#!/usr/bin/env python3
"""Bit-exact comparison of fontTools' cu2qu with kurbo's port of it.

  kdiff_py.py gen <inputs> [cubics.json]   write the fixed input set
  kdiff_py.py run <inputs> <outputs>       convert with the fontTools on
                                           PYTHONPATH, write results

Formats, with every float as the 16 hex digits of its bits:
  inputs   S <tol> <8 coordinates>            one cubic
           M <k> <tol> <8k coordinates>       k cubics converted together
  outputs  - if no approximation; else the points' coordinates, splines of a
           group separated by " | "

The Rust side (kdiff/src/main.rs) reads the same inputs and writes the same
output format.
"""
import json, math, random, struct, sys


def bits(x):
    return struct.pack(">d", x).hex()


def unbits(s):
    return struct.unpack(">d", bytes.fromhex(s))[0]


def Phi(n):
    g = (n - 2) / (4 * (n - 1))
    error = (4 / 27) * (1 + g) * (1 + g) * (1 + g)
    if n > 3:
        c = (n - 3) / (4 * (n - 1))
        w = c / 3
        error = max(error, c + 2 * w * math.sqrt(w))
    return error


def gen(path, cubics_path=None):
    rnd = random.Random(4242)
    out = open(path, "w")
    count = 0

    def single(c, tol):
        nonlocal count
        count += 1
        out.write("S %s %s\n" % (bits(float(tol)), " ".join(bits(float(v)) for p in c for v in p)))

    def multi(curves, tol):
        nonlocal count
        count += 1
        out.write("M %d %s %s\n" % (len(curves), bits(float(tol)), " ".join(bits(float(v)) for c in curves for p in c for v in p)))

    # random float cubics over scales and tolerances
    for i in range(120000):
        scale = 10 ** rnd.uniform(0, 4)
        c = [(rnd.uniform(-scale, scale), rnd.uniform(-scale, scale)) for _ in range(4)]
        single(c, rnd.choice([0.05, 0.25, 0.5, 1.0, 2.0, 2.048]))
    # random integer cubics, as in fonts
    for i in range(120000):
        size = rnd.choice([20, 60, 200, 600, 2048])
        x0, y0 = rnd.randint(-1000, 1000), rnd.randint(-1000, 1000)
        c = [(x0 + rnd.randint(-size, size), y0 + rnd.randint(-size, size)) for _ in range(4)]
        single(c, rnd.choice([0.5, 1.0, 2.048]))
    # third difference at, or a hair off, a limit of the closed form
    for i in range(40000):
        tol = rnd.choice([0.5, 1.0, 2.0])
        n = rnd.choice([2, 2, 3, 3, 4, 5, 6, 8, 12, 20, 40, 99, 100])
        limit = tol * n**3 / Phi(n)
        mag = limit * rnd.choice([1.0, 1.0, 1 - 1e-15, 1 + 1e-15, 1 - 1e-12, 1 + 1e-12, 1 - 1e-9, 1 + 1e-9, 1 - 1e-6, 1 + 1e-6, 0.999, 1.001])
        d = rnd.choice([(mag, 0.0), (0.0, mag), (-mag, 0.0), (0.0, -mag), (mag * 0.6, mag * 0.8), (-mag * 0.8, mag * 0.6)])
        span = max(50, int(limit))
        c0 = (rnd.randint(-500, 500), rnd.randint(-500, 500))
        c1 = (c0[0] + rnd.randint(-span, span), c0[1] + rnd.randint(-span, span))
        c2 = (c1[0] + rnd.randint(-span, span), c1[1] + rnd.randint(-span, span))
        c3 = (c0[0] - 3 * c1[0] + 3 * c2[0] + d[0], c0[1] - 3 * c1[1] + 3 * c2[1] + d[1])
        single([c0, c1, c2, c3], tol)
    # nearly quadratic cubics: where a single quadratic is in question
    for i in range(60000):
        size = rnd.choice([30, 100, 400, 1000])
        p0 = (rnd.randint(-size, size), rnd.randint(-size, size))
        q = (rnd.randint(-size, size), rnd.randint(-size, size))
        p3 = (rnd.randint(-size, size), rnd.randint(-size, size))
        j = rnd.choice([1, 3, 10, 30])
        c1 = (p0[0] + round(2 * (q[0] - p0[0]) / 3) + rnd.randint(-j, j), p0[1] + round(2 * (q[1] - p0[1]) / 3) + rnd.randint(-j, j))
        c2 = (p3[0] + round(2 * (q[0] - p3[0]) / 3) + rnd.randint(-j, j), p3[1] + round(2 * (q[1] - p3[1]) / 3) + rnd.randint(-j, j))
        single([p0, c1, c2, p3], rnd.choice([0.5, 1.0, 2.048]))
    # several curves at once
    for i in range(40000):
        k = rnd.randint(1, 4)
        size = rnd.choice([30, 100, 400, 2048])
        curves = [[(rnd.randint(-size, size), rnd.randint(-size, size)) for _ in range(4)] for _ in range(k)]
        multi(curves, rnd.choice([0.5, 1.0, 2.048, 0.05]))
    # cubics from fonts
    if cubics_path:
        cubics = json.load(open(cubics_path))
        for c in cubics:
            for tol in (1.0, 0.5):
                single(c, tol)
        for a, b in zip(cubics[0:40000:2], cubics[1:40000:2]):
            multi([a, b], 1.0)
    out.close()
    print("inputs:", count, file=sys.stderr)


def run(inputs, outputs):
    from fontTools.cu2qu import cu2qu
    from fontTools.cu2qu.cu2qu import curve_to_quadratic, curves_to_quadratic

    print("cu2qu from", cu2qu.__file__, "compiled:", cu2qu.COMPILED, file=sys.stderr)
    out = open(outputs, "w")

    def fmt(spline):
        return " ".join(bits(v) for p in spline for v in p)

    for line in open(inputs):
        f = line.split()
        if f[0] == "S":
            v = [unbits(s) for s in f[1:]]
            c = [(v[1], v[2]), (v[3], v[4]), (v[5], v[6]), (v[7], v[8])]
            try:
                out.write(fmt(curve_to_quadratic(c, v[0])) + "\n")
            except Exception:
                out.write("-\n")
        else:
            k = int(f[1])
            v = [unbits(s) for s in f[2:]]
            curves = [[(v[1 + 8 * i + 2 * j], v[2 + 8 * i + 2 * j]) for j in range(4)] for i in range(k)]
            try:
                out.write(" | ".join(fmt(s) for s in curves_to_quadratic(curves, [v[0]] * k)) + "\n")
            except Exception:
                out.write("-\n")
    out.close()


if __name__ == "__main__":
    if sys.argv[1] == "gen":
        gen(*sys.argv[2:])
    else:
        run(*sys.argv[2:])
