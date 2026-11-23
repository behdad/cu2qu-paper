"""Tables over the per-font results of stage 2.

usage: summarize.py RESULTS_DIR [POPULATION] [inside|all]

POPULATION: exact (default), near, flag.
inside (default): samples where no master has a negative weight; all: every
sample, including the corners of the design space that have no master.
"""

import glob
import json
import os
import sys

POPS = {
    "exact": "exactly smooth in all masters",
    "near": "smooth up to rounding",
    "flag": "marked smooth in all masters",
}


def main():
    rdir = sys.argv[1]
    pop = POPS[sys.argv[2] if len(sys.argv) > 2 else "exact"]
    which = sys.argv[3] if len(sys.argv) > 3 else "inside"
    key = "inside" if which == "inside" else "all samples"
    rows = []
    for path in sorted(glob.glob(os.path.join(rdir, "*.json"))):
        with open(path) as f:
            r = json.load(f)
        if pop not in r["summary"]:
            continue
        s = r["summary"][pop]
        if key not in s:
            if key == "all samples":
                s[key] = s["inside"]
            else:
                continue
        rows.append((os.path.basename(path)[:-5], r))
    rows.sort(key=lambda x: -x[1]["cubics"])

    print("Table 1. Joints %s; instances: %s" % (pop, key))
    print("How much the angle between the two handles of a joint changes in instances, in degrees, as counts of joints.")
    print("'c' = quadratic instance against cubic instance (what cu2qu does); 'r' = rounded quadratic against quadratic (what rounding to integers does);")
    print("'units' = the largest c as a distance: how far the nearer off-curve point is from where it would leave the angle alone; c>r = joints where cu2qu does more than rounding.")
    hdr = "%-22s %2s %8s %7s %6s | %8s %7s | %6s %6s %6s %7s %6s | %7s %7s %7s %7s | %5s" % (
        "font", "M", "cubics", "single", "%", "joints", "by-sgl",
        "c>0.1", "c>0.5", "c>1", "c max", "units",
        "r>0.1", "r>0.5", "r>1", "r max", "c>r")
    print(hdr)
    tot = dict(cubics=0, single=0, joints=0, bysgl=0, c01=0, c05=0, c1=0, r01=0, r05=0, r1=0, worse=0)
    cmax = cumax = 0.0
    for name, r in rows:
        s = r["summary"][pop]
        sm = s[key]
        single = r["n_hist"].get("1", 0)
        c = sm["change by cu2qu (all)"]
        rr = sm["change by rounding (all)"]
        cu = sm.get("change by cu2qu (next to single), units", {})
        if not c.get("count"):
            continue
        print("%-22s %2d %8d %7d %5.1f%% | %8d %7d | %6d %6d %6d %7.3f %6.3f | %7d %7d %7d %7.2f | %5d" % (
            name, r["masters"], r["cubics"], single, 100.0 * single / max(r["cubics"], 1),
            s["joints"], s["next_to_single"],
            c["gt0.1"], c["gt0.5"], c["gt1"], c["p100"], cu.get("p100", 0.0),
            rr["gt0.1"], rr["gt0.5"], rr["gt1"], rr["p100"],
            sm.get("cu2qu more than rounding", 0)))
        tot["cubics"] += r["cubics"]
        tot["single"] += single
        tot["joints"] += s["joints"]
        tot["bysgl"] += s["next_to_single"]
        tot["c01"] += c["gt0.1"]
        tot["c05"] += c["gt0.5"]
        tot["c1"] += c["gt1"]
        tot["r01"] += rr["gt0.1"]
        tot["r05"] += rr["gt0.5"]
        tot["r1"] += rr["gt1"]
        tot["worse"] += sm.get("cu2qu more than rounding", 0)
        cmax = max(cmax, c["p100"])
        cumax = max(cumax, cu.get("p100", 0.0))
    print("%-22s %2s %8d %7d %5.1f%% | %8d %7d | %6d %6d %6d %7.3f %6.3f | %7d %7d %7d %7s | %5d" % (
        "TOTAL", "", tot["cubics"], tot["single"], 100.0 * tot["single"] / max(tot["cubics"], 1),
        tot["joints"], tot["bysgl"], tot["c01"], tot["c05"], tot["c1"], cmax, cumax,
        tot["r01"], tot["r05"], tot["r1"], "", tot["worse"]))

    print()
    print("Table 2. Kinks at those joints in instances of the cubic masters, the quadratic masters and the rounded quadratic masters: joints over 0.1, 0.5 and 1 degree, and the largest")
    print("%-22s | %7s %7s %7s %7s | %7s %7s %7s %7s | %7s %7s %7s %7s" % (
        "font", "cub>0.1", ">0.5", ">1", "max", "qu>0.1", ">0.5", ">1", "max", "rnd>0.1", ">0.5", ">1", "max"))
    t = [0] * 9
    for name, r in rows:
        sm = r["summary"][pop][key]
        a, b, c = (sm[k] for k in ("cubic instances: kink", "quadratic instances: kink", "rounded quadratic: kink"))
        if not a.get("count"):
            continue
        print("%-22s | %7d %7d %7d %7.2f | %7d %7d %7d %7.2f | %7d %7d %7d %7.2f" % (
            name, a["gt0.1"], a["gt0.5"], a["gt1"], a["p100"], b["gt0.1"], b["gt0.5"], b["gt1"], b["p100"],
            c["gt0.1"], c["gt0.5"], c["gt1"], c["p100"]))
        for i, d in enumerate((a, b, c)):
            for k2, kk in enumerate(("gt0.1", "gt0.5", "gt1")):
                t[i * 3 + k2] += d[kk]
    print("%-22s | %7d %7d %7d %7s | %7d %7d %7d %7s | %7d %7d %7d %7s" % (
        "TOTAL", t[0], t[1], t[2], "", t[3], t[4], t[5], "", t[6], t[7], t[8], ""))

    print()
    print("Table 3. Candidate rules (samples inside): extra points, as %% of all points, then joints %s that cu2qu still changes by more than 0.1 degree" % pop)
    labels = None
    totals = {}
    allpoints = 0
    for name, r in rows:
        if "rules" not in r:
            continue
        rules = r["rules"]
        if labels is None:
            labels = [k for k in rules if k not in ("deviation", "linear single")]
            print("%-22s %9s | " % ("font", "points") + " | ".join("%-21s" % l[:21] for l in labels))
        cells = []
        allpoints += r["points"]
        for l in labels:
            v = rules[l]
            left = v[pop][0] if pop in v else 0
            cells.append("%7d %6.3f%% %5d" % (v["extra_points"], v["percent"], left))
            tt = totals.setdefault(l, [0, 0])
            tt[0] += v["extra_points"]
            tt[1] += left
        print("%-22s %9d | " % (name, r["points"]) + " | ".join(cells))
    if labels:
        print("%-22s %9d | " % ("TOTAL", allpoints) + " | ".join(
            "%7d %6.3f%% %5d" % (totals[l][0], 100.0 * totals[l][0] / allpoints, totals[l][1]) for l in labels))

    print()
    print("Table 4. The deviation that the rules test, for the single quadratics as given: over pairs of masters and the two ends,")
    print("the distance of the off-curve point halfway between two masters from the end tangent of the cubic halfway between them, in font units")
    print("%-22s %8s %8s | %8s %8s %8s %8s %8s" % ("font", "single", "zero", "p50", "p90", "p99", "p99.9", "max"))
    for name, r in rows:
        d = r.get("rules", {}).get("deviation")
        if not d:
            continue
        print("%-22s %8d %8d | %8.4f %8.4f %8.4f %8.4f %8.4f" % (name, d["count"], d["zero"], d["p50"], d["p90"], d["p99"], d["p99.9"], d["max"]))

    print()
    print("Table 5. Joints %s next to a single quadratic: change of the angle in the masters themselves and in instances," % pop)
    print("for the intersection of the end tangents (as is) and for a single quadratic that is linear in the control points, (3 (c1 + c2) - c0 - c3) / 4")
    print("%-22s %8s | %-27s | %-27s | %-27s | %-27s" % ("font", "joints", "as is, masters >0.1 >1 max", "as is, instances", "linear, masters", "linear, instances"))
    tot = {}
    tj = 0
    for name, r in rows:
        ls = r.get("rules", {}).get("linear single", {}).get(pop)
        if not ls:
            continue
        joints = r["summary"][pop]["next_to_single"]
        tj += joints
        cells = []
        for label in ("intersection, masters", "intersection, instances", "linear, masters", "linear, instances"):
            v = ls[label]
            cells.append("%7d %7d %9.3f  " % (v[0], v[2], v[3]))
            t = tot.setdefault(label, [0, 0, 0.0])
            t[0] += v[0]
            t[1] += v[2]
            t[2] = max(t[2], v[3])
        print("%-22s %8d | %s" % (name, joints, " | ".join(cells)))
    if tot:
        print("%-22s %8d | %s" % ("TOTAL", tj, " | ".join(
            "%7d %7d %9.3f  " % tuple(tot[l]) for l in ("intersection, masters", "intersection, instances", "linear, masters", "linear, instances"))))


if __name__ == "__main__":
    main()
