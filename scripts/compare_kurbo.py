#!/usr/bin/env python3
"""Compare kurbo's port of cu2qu with fontTools', bit for bit.

Usage: compare_kurbo.py KURBO_REPO KURBO_REV FONTTOOLS_REPO FONTTOOLS_REV cubics.json WORKDIR [TOOLCHAIN]

`CubicBez::approx_spline` of kurbo at KURBO_REV (exported and built as in
perf_kurbo.py) and `curve_to_quadratic` of fontTools at FONTTOOLS_REV (pure
Python, read out of its repository) convert every cubic of the corpus at
tolerances 1, 0.5 and 0.1. The results are compared as bits: the number of
points, and both coordinates of every point.

Then one cubic whose third difference is within rounding of a limit, where
the two may differ by a quadratic: kurbo takes the length of the third
difference as sqrt(x*x + y*y) and Python as hypot(x, y), and those can differ
in the last bit.
"""

import json
import os
import struct
import subprocess
import sys

import compare_revs
import perf_kurbo

MAIN = r"""
use std::io::{self, BufWriter, Write};
use std::{env, fs};

use kurbo::CubicBez;

fn main() {
    let mut args = env::args().skip(1);
    let path = args.next().unwrap();
    let accuracy: f64 = args.next().unwrap().parse().unwrap();
    let text = fs::read_to_string(path).unwrap();
    let out = io::stdout();
    let mut out = BufWriter::new(out.lock());
    for line in text.lines() {
        let v: Vec<f64> = line.split_whitespace().map(|x| x.parse().unwrap()).collect();
        let cubic = CubicBez::new((v[0], v[1]), (v[2], v[3]), (v[4], v[5]), (v[6], v[7]));
        match cubic.approx_spline(accuracy) {
            Some(spline) => {
                let words: Vec<String> = spline
                    .points()
                    .iter()
                    .map(|p| format!("{:016x}{:016x}", p.x.to_bits(), p.y.to_bits()))
                    .collect();
                writeln!(out, "{}", words.join(" ")).unwrap();
            }
            None => writeln!(out, "none").unwrap(),
        }
    }
}
"""

TOLERANCES = (1.0, 0.5, 0.1)
# third difference 128 to within rounding, at a tolerance of 1
AT_A_LIMIT = [[(0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (73.3702892310282, 104.88470173554994)]]


def bits(x):
    return "%016x" % struct.unpack(">Q", struct.pack(">d", x))[0]


def write_cubics(path, cubics):
    with open(path, "w") as f:
        for c in cubics:
            f.write(" ".join("%r" % float(v) for p in c for v in p) + "\n")


def kurbo(binary, path, tolerance):
    return subprocess.check_output([binary, path, repr(tolerance)], text=True).split("\n")[:-1]


def fonttools(module, cubics, tolerance):
    out = []
    for c in cubics:
        spline = module.curve_to_quadratic([tuple(p) for p in c], tolerance)
        out.append(" ".join(bits(float(x)) + bits(float(y)) for x, y in spline))
    return out


def main():
    kurbo_repo, kurbo_rev, ft_repo, ft_rev, cubics_path, work = sys.argv[1:7]
    toolchain = sys.argv[7] if len(sys.argv) > 7 else "1.99.0"
    os.makedirs(work, exist_ok=True)
    tree = os.path.join(work, "kurbo-compare")
    perf_kurbo.export(kurbo_repo, kurbo_rev, tree)
    crate = os.path.join(work, "kcompare")
    os.makedirs(os.path.join(crate, "src"), exist_ok=True)
    open(os.path.join(crate, "Cargo.toml"), "w").write(perf_kurbo.CARGO.replace("kbench", "kcompare") % os.path.join(tree, "kurbo"))
    open(os.path.join(crate, "src", "main.rs"), "w").write(MAIN)
    env = dict(os.environ, CARGO_TARGET_DIR=os.path.join(crate, "target"))
    subprocess.check_call(["cargo", "+" + toolchain, "build", "--release", "--offline", "-q"], cwd=crate, env=env)
    binary = os.path.join(crate, "target", "release", "kcompare")
    module, ft_full = compare_revs.load_revision(ft_repo, ft_rev, "cu2qu_after")
    kurbo_full = subprocess.check_output(["git", "-C", kurbo_repo, "rev-parse", kurbo_rev]).decode().strip()
    print("kurbo:     %s\nfontTools: %s (pure Python)" % (kurbo_full, ft_full))

    cubics = json.load(open(cubics_path))
    corpus = os.path.join(work, "compare-cubics.txt")
    write_cubics(corpus, cubics)
    print("%d cubics" % len(cubics))
    for tolerance in TOLERANCES:
        a = kurbo(binary, corpus, tolerance)
        b = fonttools(module, cubics, tolerance)
        assert len(a) == len(b) == len(cubics)
        count = sum(1 for x, y in zip(a, b) if len(x.split()) != len(y.split()))
        other = sum(1 for x, y in zip(a, b) if x != y) - count
        points = sum(len(x.split()) for x in a)
        print("tolerance %g: %d points; splines with a different number of points: %d; with the same number and a different bit: %d" % (tolerance, points, count, other))

    limit = os.path.join(work, "compare-limit.txt")
    write_cubics(limit, AT_A_LIMIT)
    print("\nwithin rounding of a limit, at a tolerance of 1:")
    a = kurbo(binary, limit, 1.0)
    b = fonttools(module, AT_A_LIMIT, 1.0)
    for c, x, y in zip(AT_A_LIMIT, a, b):
        dx, dy = c[3][0] - 3 * c[2][0] + 3 * c[1][0] - c[0][0], c[3][1] - 3 * c[2][1] + 3 * c[1][1] - c[0][1]
        print("  %s" % (c,))
        print("  third difference: sqrt(x*x + y*y) = %r, hypot(x, y) = %r" % ((dx * dx + dy * dy) ** 0.5, __import__("math").hypot(dx, dy)))
        print("  kurbo %d quadratics, fontTools %d" % (len(x.split()) - 2, len(y.split()) - 2))


if __name__ == "__main__":
    main()
