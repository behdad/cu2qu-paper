#!/usr/bin/env python3
"""Instructions per `CubicBez::approx_spline` call for two revisions of kurbo.

Usage: perf_kurbo.py REPO REV_BEFORE REV_AFTER cubics.json WORKDIR [TOOLCHAIN]

REPO is a kurbo git repository (nothing is checked out in it: the two
revisions are exported with `git archive` into WORKDIR). A small program is
built against each, in release mode, that converts every cubic of the corpus
at one accuracy, a given number of times. Instructions per call are the count
for 20 passes minus the count for none, over the number of calls; `perf`
counts them. Needs cargo (TOOLCHAIN, by default 1.99.0, as `cargo +TOOLCHAIN`)
and the crates kurbo depends on in the local registry.
"""

import io
import json
import os
import shutil
import subprocess
import sys
import tarfile

MAIN = r"""
use std::{env, fs, time::Instant};

use kurbo::CubicBez;

fn main() {
    let mut args = env::args().skip(1);
    let path = args.next().unwrap();
    let accuracy: f64 = args.next().unwrap().parse().unwrap();
    let passes: usize = args.next().unwrap().parse().unwrap();
    let text = fs::read_to_string(path).unwrap();
    let cubics: Vec<CubicBez> = text
        .lines()
        .map(|line| {
            let v: Vec<f64> = line.split_whitespace().map(|x| x.parse().unwrap()).collect();
            CubicBez::new((v[0], v[1]), (v[2], v[3]), (v[4], v[5]), (v[6], v[7]))
        })
        .collect();
    let mut points = 0usize;
    let mut best = f64::INFINITY;
    for _ in 0..passes {
        let start = Instant::now();
        for cubic in &cubics {
            if let Some(spline) = std::hint::black_box(cubic).approx_spline(accuracy) {
                points += spline.points().len();
            }
        }
        best = best.min(start.elapsed().as_secs_f64());
    }
    // cubics, passes, points made, nanoseconds per call in the fastest pass
    println!("{} {} {} {:.1}", cubics.len(), passes, points, best / cubics.len() as f64 * 1e9);
}
"""

CARGO = """[package]
name = "kbench"
version = "0.0.0"
edition = "2021"
publish = false

[dependencies]
kurbo = { path = "%s" }

[workspace]
"""

PASSES = 20


def export(repo, rev, dest):
    if os.path.exists(dest):
        shutil.rmtree(dest)
    os.makedirs(dest)
    data = subprocess.check_output(["git", "-C", repo, "archive", rev])
    tarfile.open(fileobj=io.BytesIO(data)).extractall(dest, filter="data")


def instructions(binary, corpus, accuracy, passes):
    best = None
    out = None
    for _ in range(3):
        p = subprocess.run(["perf", "stat", "-x,", "-e", "instructions:u", binary, corpus, str(accuracy), str(passes)], capture_output=True, text=True, check=True)
        count = int(next(line for line in p.stderr.split("\n") if "instructions" in line).split(",")[0])
        best = count if best is None else min(best, count)
        out = p.stdout.split()
    return best, out


def main():
    repo, rev_a, rev_b, cubics, work = sys.argv[1:6]
    toolchain = sys.argv[6] if len(sys.argv) > 6 else "1.99.0"
    os.makedirs(work, exist_ok=True)
    corpus = os.path.join(work, "cubics.txt")
    with open(corpus, "w") as f:
        for c in json.load(open(cubics)):
            f.write(" ".join("%r" % float(v) for p in c for v in p) + "\n")
    binaries = {}
    for label, rev in (("before", rev_a), ("after", rev_b)):
        tree = os.path.join(work, "kurbo-" + label)
        export(repo, rev, tree)
        crate = os.path.join(work, "kbench-" + label)
        os.makedirs(os.path.join(crate, "src"), exist_ok=True)
        open(os.path.join(crate, "Cargo.toml"), "w").write(CARGO % os.path.join(tree, "kurbo"))
        open(os.path.join(crate, "src", "main.rs"), "w").write(MAIN)
        env = dict(os.environ, CARGO_TARGET_DIR=os.path.join(crate, "target"))
        subprocess.check_call(["cargo", "+" + toolchain, "build", "--release", "--offline", "-q"], cwd=crate, env=env)
        binaries[label] = os.path.join(crate, "target", "release", "kbench")
        full = subprocess.check_output(["git", "-C", repo, "rev-parse", rev]).decode().strip()
        print("%s: %s" % (label, full))
    print(subprocess.check_output(["cargo", "+" + toolchain, "--version"]).decode().strip())
    print()
    print("%-8s %-10s %14s %12s %16s" % ("", "accuracy", "instructions", "ns per call", "quadratics/cubic"))
    results = {}
    for accuracy in (1.0, 0.5):
        for label in ("before", "after"):
            base, _ = instructions(binaries[label], corpus, accuracy, 0)
            count, out = instructions(binaries[label], corpus, accuracy, PASSES)
            n, passes, points, ns = int(out[0]), int(out[1]), int(out[2]), float(out[3])
            per_call = (count - base) / (n * passes)
            results[label, accuracy] = per_call
            # a spline of q quadratics has q + 2 points
            print("%-8s %-10g %14.0f %12.1f %16.4f" % (label, accuracy, per_call, ns, points / (n * passes) - 2))
    print()
    for accuracy in (1.0, 0.5):
        a, b = results["before", accuracy], results["after", accuracy]
        print("accuracy %g: %+.1f%%" % (accuracy, 100 * (b / a - 1)))


if __name__ == "__main__":
    main()
