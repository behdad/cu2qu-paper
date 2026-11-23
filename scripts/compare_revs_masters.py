#!/usr/bin/env python3
"""Compare two revisions of fontTools' cu2qu on groups of compatible masters.

Usage: compare_revs_masters.py REPO REV_BEFORE REV_AFTER MASTERS_DIR [--jobs N]

MASTERS_DIR holds the masters of some variable fonts as arrays (NAME.npz
with NAME.npz.meta.json beside it; kinks.py describes the format). Every
cubic segment is converted across its masters, with curves_to_quadratic at a
thousandth of the em, by both revisions (read out of the git repository REPO,
pure Python), and the results are compared.
"""

import glob
import json
import multiprocessing
import os
import sys

import numpy as np

import compare_revs


def work(job):
    repo, rev_a, rev_b, path, start, stop = job
    before, _ = compare_revs.load_revision(repo, rev_a, "cu2qu_before")
    after, _ = compare_revs.load_revision(repo, rev_b, "cu2qu_after")
    data = np.load(path)
    upem = json.load(open(path + ".meta.json"))["upem"]
    pts, kind = data["pts"], data["kind"]
    masters = pts.shape[1]
    tolerances = [0.001 * upem] * masters
    groups = differ = 0
    changes = {}
    for s in range(start, stop):
        if kind[s] != 1:
            continue
        curves = [[(float(z.real), float(z.imag)) for z in pts[s, m]] for m in range(masters)]
        a = before.curves_to_quadratic(curves, tolerances)
        b = after.curves_to_quadratic(curves, tolerances)
        groups += 1
        if a != b:
            differ += 1
            key = (len(a[0]) - 2, len(b[0]) - 2)
            changes[key] = changes.get(key, 0) + 1
    return os.path.basename(path)[: -len(".npz")], groups, differ, changes


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    repo, rev_a, rev_b, directory = args[:4]
    jobs = multiprocessing.cpu_count()
    if "--jobs" in sys.argv:
        jobs = int(sys.argv[sys.argv.index("--jobs") + 1])
    tasks = []
    for path in sorted(glob.glob(os.path.join(directory, "*.npz"))):
        n = len(np.load(path)["kind"])
        step = 50000
        tasks += [(repo, rev_a, rev_b, path, i, min(i + step, n)) for i in range(0, n, step)]
    totals = {}
    with multiprocessing.Pool(jobs) as pool:
        for name, groups, differ, changes in pool.imap_unordered(work, tasks):
            t = totals.setdefault(name, [0, 0, {}])
            t[0] += groups
            t[1] += differ
            for key, v in changes.items():
                t[2][key] = t[2].get(key, 0) + v
    print("before: %s\nafter:  %s" % (rev_a, rev_b))
    print("%-24s %10s %8s  %s" % ("source", "groups", "differ", "quadratics before -> after"))
    all_groups = all_differ = 0
    all_changes = {}
    for name in sorted(totals):
        groups, differ, changes = totals[name]
        all_groups += groups
        all_differ += differ
        for key, v in changes.items():
            all_changes[key] = all_changes.get(key, 0) + v
        print("%-24s %10d %8d  %s" % (name, groups, differ, ", ".join("%d -> %d: %d" % (a, b, v) for (a, b), v in sorted(changes.items()))))
    print("%-24s %10d %8d  %s" % ("total", all_groups, all_differ, ", ".join("%d -> %d: %d" % (a, b, v) for (a, b), v in sorted(all_changes.items()))))


if __name__ == "__main__":
    main()
