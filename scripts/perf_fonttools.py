#!/usr/bin/env python3
"""Instructions and time per cu2qu call for two revisions of fontTools,
compiled with Cython and in pure Python.

Usage: perf_fonttools.py REPO REV_BEFORE REV_AFTER cubics.json WORKDIR PYTHON

REPO is a fontTools git repository (nothing is checked out in it: the two
revisions are exported with `git archive` into WORKDIR). PYTHON is an
interpreter that has Cython and fontTools' requirements installed. Needs
`perf` (user-space instruction counts).

For each tree and workload (see cu2qu_workload.py) this prints instructions
per call, the workload's count minus the load-only baseline, least of three
runs; and microseconds per call, best of seven.
"""

import io
import os
import shutil
import subprocess
import sys
import tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
WORKLOAD = os.path.join(HERE, "cu2qu_workload.py")
MODES = ("single", "half", "multi")


def export(repo, rev, dest):
    if os.path.exists(dest):
        shutil.rmtree(dest)
    os.makedirs(dest)
    data = subprocess.check_output(["git", "-C", repo, "archive", rev])
    tarfile.open(fileobj=io.BytesIO(data)).extractall(dest, filter="data")


def build(tree, python):
    env = dict(os.environ, FONTTOOLS_WITH_CYTHON="1")
    log = tree + ".build.log"
    with open(log, "w") as f:
        subprocess.check_call([python, "setup.py", "-q", "build_ext", "--inplace"], cwd=tree, env=env, stdout=f, stderr=f)


def run(python, lib, cubics, mode, *extra):
    env = dict(os.environ, PYTHONPATH=lib, PYTHONHASHSEED="0")
    return subprocess.run([python, WORKLOAD, cubics, mode, *extra], env=env, cwd="/", capture_output=True, text=True, check=True).stdout.split("\n")


def instructions(python, lib, cubics, mode):
    env = dict(os.environ, PYTHONPATH=lib, PYTHONHASHSEED="0")
    best = None
    calls = compiled = None
    for _ in range(3):
        p = subprocess.run(
            ["perf", "stat", "-x,", "-e", "instructions:u", python, WORKLOAD, cubics, mode],
            env=env, cwd="/", capture_output=True, text=True, check=True,
        )
        count = int(next(line for line in p.stderr.split("\n") if "instructions" in line).split(",")[0])
        best = count if best is None else min(best, count)
        _, calls, compiled = p.stdout.strip().split("\n")[-1].split()
    return best, int(calls), compiled


def main():
    repo, rev_a, rev_b, cubics, work, python = sys.argv[1:7]
    cubics = os.path.abspath(cubics)
    os.makedirs(work, exist_ok=True)
    trees = []
    for label, rev in (("before", rev_a), ("after", rev_b)):
        full = subprocess.check_output(["git", "-C", repo, "rev-parse", rev]).decode().strip()
        pure = os.path.join(work, "ft-%s-pure" % label)
        cy = os.path.join(work, "ft-%s-cy" % label)
        export(repo, rev, pure)
        export(repo, rev, cy)
        build(cy, python)
        trees.append((label + ", compiled", cy, full))
        trees.append((label + ", pure Python", pure, full))
        print("%s: %s" % (label, full))
    print(subprocess.check_output([python, "-c", "import sys, Cython; print('Python', sys.version.split()[0], 'Cython', Cython.__version__)"]).decode().strip())
    print()
    print("%-22s %-9s %14s %14s %14s" % ("", "", *MODES))
    results = {}
    for label, tree, _ in trees:
        lib = os.path.join(tree, "Lib")
        base, _, compiled = instructions(python, lib, cubics, "none")
        want = "True" if "compiled" in label else "False"
        assert compiled == want, (label, compiled)
        row = []
        times = []
        for mode in MODES:
            count, calls, _ = instructions(python, lib, cubics, mode)
            row.append((count - base) / calls)
            times.append(float(run(python, lib, cubics, mode, "time")[0]))
        results[label] = (row, times)
        print("%-22s %-9s %14.0f %14.0f %14.0f" % (label, "instr.", *row))
        print("%-22s %-9s %14.2f %14.2f %14.2f" % (label, "microsec.", *times))
    print()
    for kind in ("compiled", "pure Python"):
        a = results["before, " + kind]
        b = results["after, " + kind]
        print("%-22s %-9s %14s %14s %14s" % ("change, " + kind, "instr.", *("%+.1f%%" % (100 * (y / x - 1)) for x, y in zip(a[0], b[0]))))
        print("%-22s %-9s %14s %14s %14s" % ("change, " + kind, "microsec.", *("%+.1f%%" % (100 * (y / x - 1)) for x, y in zip(a[1], b[1]))))


if __name__ == "__main__":
    main()
