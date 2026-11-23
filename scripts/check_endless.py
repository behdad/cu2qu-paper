#!/usr/bin/env python3
"""Inputs on which cu2qu's fit test does not return.

Usage: check_endless.py REPO REV_BEFORE REV_AFTER

Each input is converted in a process of its own, with the two revisions of
fontTools' cu2qu read out of the git repository REPO (pure Python), and with
whatever fontTools the interpreter has installed (compiled, if it came from a
wheel). The outcome is the result, the exception, or the signal that killed
the process.
"""

import os
import signal
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

INPUTS = [
    ("[(0, 0), (0, 0), (0, 0), (25, 0)]", "100 / 27"),
    ("[(-45, -3), (-30, -2), (0, 5), (0, 9)]", "20 / 3"),
    # in the search over n, at a limit of the spline: |d3| = 128 x 5, and 843.75 x 4
    ("[(-258, -394), (-148, -194), (-155, -13), (105, 661)]", "5"),
    ("[(-22, 289), (31, 324), (-137, 476), (1499, 3445)]", "4"),
    ("[(0, 0), (1, 1), (float('nan'), 0), (3, 0)]", "1.0"),
    ("[(0, 0), (30, 40), (70, 40), (100, 0)]", "float('nan')"),
]

CHILD = r"""
import sys
sys.setrecursionlimit(10000)
if len(sys.argv) > 3:
    sys.path.insert(0, %(here)r)
    from compare_revs import load_revision
    module, _ = load_revision(sys.argv[3], sys.argv[4], "cu2qu_rev")
else:
    from fontTools.cu2qu import cu2qu as module
try:
    result = module.curve_to_quadratic(eval(sys.argv[1]), eval(sys.argv[2]))
    print("returns %%d points%%s" %% (len(result), "" if len(result) > 3 else ": " + str(result)))
except RecursionError:
    print("RecursionError")
except Exception as e:
    print(type(e).__name__)
""" % {"here": HERE}


def outcome(args):
    p = subprocess.run([sys.executable, "-c", CHILD, *args], capture_output=True, text=True)
    if p.returncode < 0:
        return "killed by %s" % signal.Signals(-p.returncode).name
    return p.stdout.strip() or ("exit %d: %s" % (p.returncode, p.stderr.strip().split("\n")[-1]))


def main():
    repo, rev_a, rev_b = sys.argv[1:4]
    import fontTools
    from fontTools.cu2qu import cu2qu

    columns = [
        ("installed fontTools %s (%s)" % (fontTools.version, "compiled" if cu2qu.COMPILED else "pure Python"), []),
        ("before %s (pure Python)" % rev_a, [repo, rev_a]),
        ("after %s (pure Python)" % rev_b, [repo, rev_b]),
    ]
    for curve, tolerance in INPUTS:
        print("curve_to_quadratic(%s, %s)" % (curve, tolerance))
        for name, extra in columns:
            print("    %-46s %s" % (name, outcome([curve, tolerance, *extra])))


if __name__ == "__main__":
    main()
