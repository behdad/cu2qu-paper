#!/usr/bin/env python3
"""Inputs on which kurbo's port of cu2qu does not return.

Usage: check_endless_kurbo.py REPO REV_BEFORE REV_AFTER WORKDIR [TOOLCHAIN]

Like check_endless.py, for `CubicBez::approx_spline` in two revisions of
kurbo (exported from the git repository REPO into WORKDIR and built in
release mode; see perf_kurbo.py). Each input is converted in a process of its
own; the outcome is the number of points, `None`, or the signal that ended
the process. A stack that overflows ends a Rust program with SIGABRT.
"""

import os
import signal
import subprocess
import sys

import perf_kurbo

MAIN = r"""
use std::env;

use kurbo::CubicBez;

fn main() {
    let v: Vec<f64> = env::args().skip(1).map(|x| x.parse().unwrap()).collect();
    let cubic = CubicBez::new((v[0], v[1]), (v[2], v[3]), (v[4], v[5]), (v[6], v[7]));
    match cubic.approx_spline(v[8]) {
        Some(spline) => println!("returns {} points", spline.points().len()),
        None => println!("returns None"),
    }
}
"""

# (cubic, accuracy, what it is)
INPUTS = [
    # the two that fontTools did not return from (check_endless.py)
    ([(0, 0), (0, 0), (0, 0), (25, 0)], "100/27", "single quadratic, fontTools' tie"),
    ([(-45, -3), (-30, -2), (0, 5), (0, 9)], "20/3", "single quadratic, fontTools' tie"),
    # the single quadratic as far as the accuracy to the last bit, in kurbo's arithmetic
    ([(-31, 97), (-51, 93), (-91, 76), (-139, 52)], "1.987615979999813", "single quadratic, kurbo's tie"),
    ([(-20, 7), (-20, -27.5), (-12, -32.5), (4, -32)], "3.5555555555555554", "single quadratic, kurbo's tie"),
    # in the search over n: third difference 843.75 x 4, the limit of six quadratics
    ([(386, -158), (527, -28), (480, -99), (2270, 2329)], "4", "spline of six at its limit"),
    # the two that fontTools fails on in the search over n (check_endless.py)
    ([(-258, -394), (-148, -194), (-155, -13), (105, 661)], "5", "spline of three at its limit"),
    ([(-22, 289), (31, 324), (-137, 476), (1499, 3445)], "4", "spline of six at its limit"),
    ([(0, 0), (float("nan"), 0), (1, 1), (2, 2)], "1", "NaN coordinate"),
    ([(0, 0), (0, 100), (100, 200), (354, 300)], "nan", "NaN accuracy"),
]


def number(text):
    if "/" in text:
        a, b = text.split("/")
        return float(a) / float(b)
    return float(text)


def outcome(binary, cubic, accuracy):
    args = [repr(float(v)) for p in cubic for v in p] + [repr(number(accuracy))]
    try:
        p = subprocess.run([binary, *args], capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return "no result in two minutes"
    if p.returncode < 0:
        what = "stack overflow, " if "overflowed its stack" in p.stderr else ""
        return "%skilled by %s" % (what, signal.Signals(-p.returncode).name)
    return p.stdout.strip() or "exit %d" % p.returncode


def main():
    repo, rev_a, rev_b, work = sys.argv[1:5]
    toolchain = sys.argv[5] if len(sys.argv) > 5 else "1.99.0"
    os.makedirs(work, exist_ok=True)
    binaries = []
    for label, rev in (("before", rev_a), ("after", rev_b)):
        tree = os.path.join(work, "kurbo-" + label)
        perf_kurbo.export(repo, rev, tree)
        crate = os.path.join(work, "kendless-" + label)
        os.makedirs(os.path.join(crate, "src"), exist_ok=True)
        open(os.path.join(crate, "Cargo.toml"), "w").write(perf_kurbo.CARGO.replace("kbench", "kendless") % os.path.join(tree, "kurbo"))
        open(os.path.join(crate, "src", "main.rs"), "w").write(MAIN)
        env = dict(os.environ, CARGO_TARGET_DIR=os.path.join(crate, "target"))
        subprocess.check_call(["cargo", "+" + toolchain, "build", "--release", "--offline", "-q"], cwd=crate, env=env)
        full = subprocess.check_output(["git", "-C", repo, "rev-parse", "--short=9", rev]).decode().strip()
        binaries.append(("%s %s" % (label, full), os.path.join(crate, "target", "release", "kendless")))
    for cubic, accuracy, what in INPUTS:
        points = ", ".join("(%s, %s)" % tuple("%g" % v for v in p) for p in cubic)
        print("CubicBez(%s).approx_spline(%s)    [%s]" % (points, accuracy, what))
        for name, binary in binaries:
            print("    %-24s %s" % (name, outcome(binary, cubic, accuracy)))


if __name__ == "__main__":
    main()
