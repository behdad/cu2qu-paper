#!/usr/bin/env python3
"""Compare two revisions of fontTools' cu2qu on a corpus of cubics, and both
with the closed form.

Usage: compare_revs.py REPO REV_BEFORE REV_AFTER cubics.json

The two revisions of Lib/fontTools/cu2qu/{cu2qu,errors}.py are read out of
the git repository REPO (nothing is checked out) and run in pure Python.
For tolerances 1, 0.5 and 0.1 this prints the quadratics per cubic of each,
the cubics whose results differ, whether the number of quadratics is the
one cu2qu_math.count() gives, and whether a single quadratic is returned
exactly when the closed form of its error (measure.py) is within tolerance.
"""

import importlib
import json
import os
import subprocess
import sys
import tempfile

import cu2qu_math as M
import measure

TOLERANCES = (1.0, 0.5, 0.1)


def load_revision(repo, rev, name):
    directory = tempfile.mkdtemp(prefix="cu2qu-rev-")
    package = os.path.join(directory, name)
    os.mkdir(package)
    open(os.path.join(package, "__init__.py"), "w").close()
    for module in ("cu2qu", "errors"):
        source = subprocess.check_output(["git", "-C", repo, "show", "%s:Lib/fontTools/cu2qu/%s.py" % (rev, module)])
        with open(os.path.join(package, module + ".py"), "wb") as f:
            f.write(source)
    sys.path.insert(0, directory)
    module = importlib.import_module(name + ".cu2qu")
    sys.path.pop(0)
    full = subprocess.check_output(["git", "-C", repo, "rev-parse", rev]).decode().strip()
    return module, full


def main():
    repo, rev_a, rev_b, path = sys.argv[1:5]
    before, full_a = load_revision(repo, rev_a, "cu2qu_before")
    after, full_b = load_revision(repo, rev_b, "cu2qu_after")
    print("before: %s\nafter:  %s" % (full_a, full_b))
    print("compiled: %s %s" % (before.COMPILED, after.COMPILED))
    cubics = [[tuple(p) for p in c] for c in json.load(open(path))]
    print("%d cubics" % len(cubics))
    for tol in TOLERANCES:
        total_a = total_b = differ = off_formula_a = off_formula_b = off_single = singles = 0
        print("\ntolerance %g" % tol)
        for c in cubics:
            a = before.curve_to_quadratic(c, tol)
            b = after.curve_to_quadratic(c, tol)
            na, nb = len(a) - 2, len(b) - 2
            total_a += na
            total_b += nb
            z = [complex(*p) for p in c]
            d3 = abs(M.third_difference(z))
            want = M.count(d3, tol)
            if na > 1 and na != want:
                off_formula_a += 1
            if nb > 1 and nb != want:
                off_formula_b += 1
            # the single quadratic: the closed form of its error against the test
            error = measure.single_equal_parameter_error(z)
            fits = error is not None and error <= tol
            singles += nb == 1
            if fits != (nb == 1):
                off_single += 1
            if a != b:
                differ += 1
                print("  differs: %s  |d3| / tolerance = %.17g: %d -> %d quadratics" % (c, d3 / tol, na, nb))
        n = len(cubics)
        print("  quadratics per cubic: before %.4f, after %.4f" % (total_a / n, total_b / n))
        print("  results that differ: %d of %d" % (differ, n))
        print("  splines whose number is not the closed form's: before %d, after %d" % (off_formula_a, off_formula_b))
        print("  single quadratics after: %d; cubics that the closed form of its error decides otherwise: %d" % (singles, off_single))


if __name__ == "__main__":
    main()
