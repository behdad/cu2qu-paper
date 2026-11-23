#!/usr/bin/env python3
"""Convert the inputs of the kurbo/fontTools comparison (the format of
kdiff_py.py: "S <tol> <8 coords>" or "M <k> <tol> <8k coords>", every float
as the 16 hex digits of its bits) with the fontTools on PYTHONPATH, and
write the results bit for bit, one line per input:

    -        no approximation (ApproxNotFoundError, or any other exception)
    R        RecursionError (pure Python: the fit test did not end)
    CRASH    the interpreter died (compiled: the fit test did not end)
    else the points' coordinates as hex bits, splines of a group separated
    by " | "

Usage: ft_run.py <inputs> <outputs>

The conversion runs in a child process that is restarted after the input
that killed it, so a crash is an outcome and not the end of the run.
"""
import os, struct, subprocess, sys


def bits(x):
    return struct.pack(">d", x).hex()


def unbits(s):
    return struct.unpack(">d", bytes.fromhex(s))[0]


def worker(inputs, outputs, start):
    from fontTools.cu2qu.cu2qu import curve_to_quadratic, curves_to_quadratic
    import resource

    # no core file when the compiled module overflows the stack
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    sys.setrecursionlimit(400)
    out = open(outputs, "a")

    def fmt(spline):
        return " ".join(bits(v) for p in spline for v in p)

    with open(inputs) as f:
        for i, line in enumerate(f):
            if i < start:
                continue
            fields = line.split()
            try:
                if fields[0] == "S":
                    v = [unbits(s) for s in fields[1:]]
                    c = [(v[1], v[2]), (v[3], v[4]), (v[5], v[6]), (v[7], v[8])]
                    r = fmt(curve_to_quadratic(c, v[0]))
                else:
                    k = int(fields[1])
                    v = [unbits(s) for s in fields[2:]]
                    curves = [[(v[1 + 8 * i2 + 2 * j], v[2 + 8 * i2 + 2 * j]) for j in range(4)] for i2 in range(k)]
                    r = " | ".join(fmt(s) for s in curves_to_quadratic(curves, [v[0]] * k))
            except RecursionError:
                r = "R"
            except Exception:
                r = "-"
            out.write(r + "\n")
            out.flush()
    out.close()


def main(inputs, outputs):
    from fontTools.cu2qu import cu2qu

    print("cu2qu from", cu2qu.__file__, "compiled:", cu2qu.COMPILED, file=sys.stderr)
    open(outputs, "w").close()
    total = sum(1 for _ in open(inputs))
    crashes = 0
    while True:
        done = sum(1 for _ in open(outputs))
        if done >= total:
            break
        r = subprocess.run([sys.executable, __file__, "--worker", inputs, outputs, str(done)], stderr=subprocess.DEVNULL)
        if r.returncode != 0:
            # the line after the last one written killed the worker
            with open(outputs, "rb+") as f:
                data = f.read()
                # drop a partial last line, if any
                if data and not data.endswith(b"\n"):
                    f.truncate(data.rfind(b"\n") + 1)
            with open(outputs, "a") as f:
                f.write("CRASH\n")
            crashes += 1
    print("inputs:", total, "crashes:", crashes, file=sys.stderr)


if __name__ == "__main__":
    if sys.argv[1] == "--worker":
        worker(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    else:
        main(sys.argv[1], sys.argv[2])
