"""Stage 1: extract compatible cubic masters and convert them with cu2qu.

usage: stage1.py NAME KIND PATH FONTNUMBER OUT.npz [JOBS] [TOL_PER_EM]

KIND is cff2 (a CFF2 variable font; masters at the region peaks), ds (a
designspace of UFOs, or a .glyphs file), ttvf (a TrueType variable font, its
quadratics raised to cubics and rounded to integers: masters at the ends of
the axes) or ttvf-exact (the same without the rounding).  The tolerance is
units-per-em times TOL_PER_EM (default 0.001, what ufo2ft and fontc use).
"""

import json
import multiprocessing
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kinklib


def work(job):
    kind, path, num, lo, hi, tol_per_em = job
    if kind == "cff2":
        arrays, locs, upem, skipped = kinklib.extract_cff2(path, num, (lo, hi))
    elif kind in ("ttvf", "ttvf-exact"):
        arrays, locs, upem, skipped = kinklib.extract_ttvf(
            path, num, (lo, hi), rounded=kind == "ttvf"
        )
    else:
        arrays, locs, upem, skipped = kinklib.extract_designspace(path, (lo, hi))
    kinklib.convert(arrays, upem * tol_per_em)
    return arrays, locs, upem, skipped


def count_glyphs(kind, path, num):
    if kind in ("cff2", "ttvf", "ttvf-exact"):
        from fontTools.ttLib import TTFont

        font = TTFont(path, fontNumber=num) if num >= 0 else TTFont(path)
        return len(font.getGlyphOrder())
    return None


def main():
    name, kind, path, num, out = sys.argv[1:6]
    num = int(num)
    jobs = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    tol_per_em = float(sys.argv[7]) if len(sys.argv) > 7 else 0.001
    t0 = time.time()
    total = count_glyphs(kind, path, num)
    if total is None or jobs == 1:
        results = [work((kind, path, num, 0, 10**9, tol_per_em))]
    else:
        step = (total + jobs - 1) // jobs
        todo = [
            (kind, path, num, lo, min(lo + step, total), tol_per_em)
            for lo in range(0, total, step)
        ]
        with multiprocessing.Pool(jobs) as pool:
            results = pool.map(work, todo)
    # Concatenate, shifting the indices.
    parts = {}
    soff = goff = 0
    skipped = {}
    for arrays, locs, upem, skip in results:
        S = len(arrays["kind"])
        prev = arrays["prev"].copy()
        prev[prev >= 0] += soff
        arrays["prev"] = prev
        arrays["gid"] = arrays["gid"] + goff
        soff += S
        goff += len(arrays["names"])
        for k, v in arrays.items():
            parts.setdefault(k, []).append(v)
        for k, v in skip.items():
            skipped[k] = skipped.get(k, 0) + v
    arrays = {k: np.concatenate(v) for k, v in parts.items()}
    meta = dict(
        name=name,
        path=path,
        font_number=num,
        locations=locs,
        upem=upem,
        tol=upem * tol_per_em,
        skipped=skipped,
        glyphs=int(len(arrays["names"])),
        segments=int(len(arrays["kind"])),
        seconds=time.time() - t0,
    )
    names = arrays.pop("names")
    np.savez_compressed(out, **arrays)
    with open(out + ".names.json", "w") as f:
        json.dump([str(x) for x in names], f)
    with open(out + ".meta.json", "w") as f:
        json.dump(meta, f, indent=1)
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
