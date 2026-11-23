"""Stage 2: measure kinks at joints in instances interpolated from masters.

usage: stage2.py DATA.npz [OUT.json] [--rules]

For every joint (a point where two segments of a contour meet) whose two
handles are not of zero length in any master, the angle between the incoming
and the outgoing handle is computed in instances interpolated

  c  from the cubic masters,
  q  from the quadratic masters as curves_to_quadratic gives them,
  l  from quadratic masters where no curve got a single quadratic (the single
     ones get the spline of two instead),
  r  from the quadratic masters (q) with every point rounded to integers.

A handle is the vector between the on-curve point and the neighbouring
control point; for a line it is the line.

Instances are sampled between the masters: on one axis at tenths of the way
between neighbouring masters, on several axes at the grid of the masters'
coordinates and the midpoints between them.  A sample where every master has a
weight that is not negative is "inside"; the others (corners of the design
space that have no master) are extrapolations and are reported apart.

With --rules, candidate rules for when curves converted together may get
single quadratics are evaluated: what they cost in points and what kinks they
leave.
"""

import itertools
import json
import os
import sys

import numpy as np

DEG = 180.0 / np.pi


def master_weights(locations):
    """Sample locations between masters and the weights of the masters there.

    Returns (list of location dicts, W[L, M])."""
    from fontTools.varLib.models import VariationModel

    axes = sorted({k for loc in locations for k in loc})
    full = [{a: loc.get(a, 0.0) for a in axes} for loc in locations]
    model = VariationModel(full, axisOrder=axes)
    per_axis = []
    for a in axes:
        vals = sorted({loc[a] for loc in full})
        cand = set(vals)
        for lo, hi in zip(vals, vals[1:]):
            if len(axes) == 1:
                fracs = [i / 10.0 for i in range(1, 10)]
            else:
                fracs = (0.5,)
            for f in fracs:
                cand.add(lo + (hi - lo) * f)
        per_axis.append(sorted(cand))
    masters = {tuple(loc[a] for a in axes) for loc in full}
    samples = []
    W = []
    for combo in itertools.product(*per_axis):
        if combo in masters:
            continue
        loc = dict(zip(axes, combo))
        samples.append(loc)
        W.append(model.getMasterScalars(loc))
    return samples, np.array(W, dtype=np.float64)


def kink(inh, outh):
    """Signed angle in degrees from the incoming to the outgoing handle."""
    return np.angle(outh * np.conj(inh)) * DEG


def wrap(a):
    return (a + 180.0) % 360.0 - 180.0


def oround(z):
    return np.floor(z.real + 0.5) + 1j * np.floor(z.imag + 0.5)


def cross(a, b):
    return a.real * b.imag - a.imag * b.real


def dot(a, b):
    return a.real * b.real + a.imag * b.imag


def load(path):
    data = dict(np.load(path))
    with open(path + ".meta.json") as f:
        meta = json.load(f)
    with open(path + ".names.json") as f:
        names = json.load(f)
    return data, meta, names


def single_fractions(pts):
    """For cubics pts[S, M, 4]: s1, s2 with q = c0 + s1 (c1 - c0)
    = c3 + s2 (c2 - c3), the intersection of the end tangents."""
    c0, c1, c2, c3 = (pts[..., k] for k in range(4))
    h1 = c1 - c0
    h2 = c2 - c3
    den = cross(h1, h2)
    with np.errstate(divide="ignore", invalid="ignore"):
        s1 = cross(c3 - c0, h2) / den
        s2 = cross(c3 - c0, h1) / den
    return s1, s2, h1, h2


def pair_deviation(s, h):
    """Largest, over pairs of masters, distance of the midpoint of the two
    masters' off-curve points from the tangent line of the midpoint cubic:
    |s_i - s_j| |h_i x h_j| / (2 |h_i + h_j|).

    s[S, M] handle factors (off-curve = on-curve + s * h), h[S, M] handles."""
    S, M = s.shape
    out = np.zeros(S)
    for i in range(M):
        for j in range(i + 1, M):
            num = np.abs(s[:, i] - s[:, j]) * np.abs(cross(h[:, i], h[:, j]))
            den = 2.0 * np.abs(h[:, i] + h[:, j])
            with np.errstate(divide="ignore", invalid="ignore"):
                d = np.where(num == 0.0, 0.0, num / den)
            out = np.fmax(out, d)
    return out


def pair_fraction_gap(s, h):
    """The same without regard to direction: how far apart, along the handle,
    the off-curve points of two masters are for handles of the same direction:
    |s_i - s_j| |h_i| |h_j| / (|h_i| + |h_j|), the largest over pairs."""
    S, M = s.shape
    out = np.zeros(S)
    for i in range(M):
        for j in range(i + 1, M):
            a = np.abs(h[:, i])
            b = np.abs(h[:, j])
            with np.errstate(divide="ignore", invalid="ignore"):
                d = np.where(a + b > 0, np.abs(s[:, i] - s[:, j]) * a * b / (a + b), 0.0)
            out = np.fmax(out, d)
    return out


def single_deviation(data, measure=pair_deviation):
    """Per segment: for the curves that got single quadratics, the largest
    `measure` of the two ends; nan elsewhere."""
    pts = data["pts"]
    single = (data["kind"] == 1) & (data["n"] == 1)
    idx = np.nonzero(single)[0]
    dev = np.full(len(single), np.nan)
    if len(idx):
        p = pts[idx]
        # From the points cu2qu gave: q = c0 + s1 h1 = c3 + s2 h2.
        h1 = p[:, :, 1] - p[:, :, 0]
        h2 = p[:, :, 2] - p[:, :, 3]
        g1 = data["qfirst"][idx] - p[:, :, 0]
        g2 = data["qlast"][idx] - p[:, :, 3]
        with np.errstate(divide="ignore", invalid="ignore"):
            s1 = np.where(h1 != 0, (g1 / h1).real, 0.0)
            s2 = np.where(h2 != 0, (g2 / h2).real, 0.0)
        dev[idx] = np.fmax(measure(s1, h1), measure(s2, h2))
    return dev


def analyze(data, meta, bump=None, only=None, midpoint=False):
    """Per joint measures.

    bump: boolean per segment, the single-quadratic curves to give two
    quadratics instead (a candidate rule); None for none.
    only: boolean per segment; if given, only the joints at the start of
    these segments are computed.
    midpoint: the single quadratics get the off-curve point that is linear in
    the control points, (3 (c1 + c2) - c0 - c3) / 4, instead of the
    intersection of the end tangents."""
    pts = data["pts"]
    kind = data["kind"]
    prev = data["prev"]
    n = data["n"].astype(np.int64)
    S, M, _ = pts.shape

    sel = prev >= 0
    if only is not None:
        sel &= only
    j = np.nonzero(sel)[0]
    p = prev[j]

    def handles(seg, at_start):
        """Handles of the segments seg at their start (out) or end (in), for
        the four variants."""
        P = pts[seg]
        cub = (kind[seg] == 1)[:, None]
        nn = n[seg]
        line = P[:, :, 1] - P[:, :, 0]
        if at_start:
            on = P[:, :, 0]
            hc = line  # c1 - c0 for a cubic as well
            with np.errstate(invalid="ignore"):
                hq = np.where(cub, data["qfirst"][seg] - on, line)
        else:
            on = np.where(cub, P[:, :, 3], P[:, :, 1])
            hc = np.where(cub, P[:, :, 3] - P[:, :, 2], line)
            with np.errstate(invalid="ignore"):
                hq = np.where(cub, on - data["qlast"][seg], line)
        hl = np.where(cub, hc * (1.5 / np.maximum(nn, 2))[:, None], hc)
        if bump is not None:
            hq = np.where(bump[seg][:, None], hl, hq)
        if midpoint:
            with np.errstate(invalid="ignore"):
                d3 = P[:, :, 3] - 3 * P[:, :, 2] + 3 * P[:, :, 1] - P[:, :, 0]
                hm = 1.5 * hc - 0.25 * d3
            hq = np.where(cub & (nn == 1)[:, None], hm, hq)
        # Rounded: the off-curve point next to the joint, and the joint; for a
        # line its two ends.
        ron = oround(on)
        if at_start:
            other = np.where(cub, on + hq, P[:, :, 1])
            hr = oround(other) - ron
        else:
            other = np.where(cub, on - hq, P[:, :, 0])
            hr = ron - oround(other)
        return dict(c=hc, q=hq, l=hl, r=hr)

    I = handles(p, False)
    O = handles(j, True)
    valid = np.all((I["c"] != 0) & (O["c"] != 0), axis=1)

    # Smoothness of the cubic masters.
    cr = cross(I["c"], O["c"])
    dt = dot(I["c"], O["c"])
    exact = np.all((cr == 0) & (dt > 0), axis=1)
    near = np.all(
        (dt > 0) & (np.abs(cr) <= 0.75 * (np.abs(I["c"]) + np.abs(O["c"]))), axis=1
    )
    with np.errstate(invalid="ignore"):
        km = np.max(np.abs(kink(I["c"], O["c"])), axis=1)
        # The masters themselves, converted: how much the angle changed.
        dqm = np.max(
            np.abs(wrap(kink(I["q"], O["q"]) - kink(I["c"], O["c"]))), axis=1
        )

    # How far the tangent turns between masters.
    theta = np.zeros(len(j))
    with np.errstate(invalid="ignore"):
        for a in range(M):
            for b in range(a + 1, M):
                t = np.abs(np.angle(O["c"][:, a] * np.conj(O["c"][:, b]))) * DEG
                theta = np.fmax(theta, t)

    samples, W = master_weights(meta["locations"])
    inside = np.min(W, axis=1) >= -1e-9
    res = {}
    with np.errstate(invalid="ignore", divide="ignore"):
        HI = {k: I[k] @ W.T for k in "cqlr"}
        HO = {k: O[k] @ W.T for k in "cqlr"}
        K = {k: kink(HI[k], HO[k]) for k in "cqlr"}
        # An instance where a handle (nearly) vanishes has no tangent there.
        scale = np.max(np.abs(I["c"]), axis=1) + np.max(np.abs(O["c"]), axis=1)
        gone = np.zeros(K["c"].shape, dtype=bool)
        for k in "cqlr":
            gone |= (np.abs(HI[k]) < 1e-6 * scale[:, None]) | (
                np.abs(HO[k]) < 1e-6 * scale[:, None]
            )
    for k in "cqlr":
        K[k] = np.where(gone, 0.0, K[k])
    dq = np.abs(wrap(K["q"] - K["c"]))
    dl = np.abs(wrap(K["l"] - K["c"]))
    dr = np.abs(wrap(K["r"] - K["q"]))
    # The same changes as a distance: how far the end of the shorter handle
    # moves sideways.
    short = np.minimum(np.abs(HI["q"]), np.abs(HO["q"]))
    short = np.where(gone, 0.0, short)
    xq = np.abs(np.sin(dq / DEG)) * short
    xr = np.abs(np.sin(dr / DEG)) * short

    def reduce(x, cols):
        if not np.any(cols):
            return np.zeros(x.shape[0])
        return np.max(x[:, cols], axis=1)

    for suffix, cols in (("", inside), ("_all", np.ones(len(inside), dtype=bool))):
        res["kc" + suffix] = reduce(np.abs(K["c"]), cols)
        res["kq" + suffix] = reduce(np.abs(K["q"]), cols)
        res["kr" + suffix] = reduce(np.abs(K["r"]), cols)
        res["dq" + suffix] = reduce(dq, cols)
        res["dl" + suffix] = reduce(dl, cols)
        res["dr" + suffix] = reduce(dr, cols)
        res["xq" + suffix] = reduce(xq, cols)
        res["xr" + suffix] = reduce(xr, cols)
    res["dq_loc"] = np.argmax(np.where(inside[None, :], dq, -1.0), axis=1)
    cubp = kind[p] == 1
    cubj = kind[j] == 1
    res.update(
        j=j,
        p=p,
        valid=valid,
        exact=exact,
        near=near,
        km=km,
        dqm=dqm,
        theta=theta,
        tin=np.where(cubp, n[p], 0),
        tout=np.where(cubj, n[j], 0),
        flag=data["flag"][j],
        samples=samples,
        inside=inside,
        W=W,
    )
    if bump is not None:
        res["tin"] = np.where(bump[p], 2, res["tin"])
        res["tout"] = np.where(bump[j], 2, res["tout"])
    return res


QS = (50, 90, 99, 99.9, 100)
THRESHOLDS = (0.1, 0.5, 1.0)


def describe(x):
    if len(x) == 0:
        return dict(count=0)
    out = dict(count=int(len(x)))
    for q, v in zip(QS, np.percentile(x, QS)):
        out["p%g" % q] = float(v)
    for t in THRESHOLDS:
        out["gt%g" % t] = int(np.sum(x > t))
    return out


def fmt(d):
    if not d.get("count"):
        return "%9d" % 0
    return "%9d  %8.4f %8.4f %8.4f %8.4f %8.4f  %8d %8d %8d" % (
        d["count"],
        d["p50"],
        d["p90"],
        d["p99"],
        d["p99.9"],
        d["p100"],
        d["gt0.1"],
        d["gt0.5"],
        d["gt1"],
    )


HEADER = "%-36s %9s  %8s %8s %8s %8s %8s  %8s %8s %8s" % (
    "",
    "joints",
    "p50",
    "p90",
    "p99",
    "p99.9",
    "max",
    ">0.1",
    ">0.5",
    ">1",
)


def populations(res):
    pops = [
        ("exactly smooth in all masters", res["valid"] & res["exact"]),
        ("smooth up to rounding", res["valid"] & res["near"]),
    ]
    if np.any(res["flag"] >= 0):
        pops.append(("marked smooth in all masters", res["valid"] & (res["flag"] == 1)))
    return pops


def report(res, title, out=sys.stdout):
    summary = {}
    single = (res["tin"] == 1) | (res["tout"] == 1)
    both_single = (res["tin"] == 1) & (res["tout"] == 1)
    one_single = single & ~both_single
    other = np.where(res["tin"] == 1, res["tout"], res["tin"])
    other_spline = other >= 2
    other_line = other == 0
    print("== %s" % title, file=out)
    print(
        "joints %d, with handles in all masters %d; samples %d, of them inside the masters %d"
        % (
            len(res["j"]),
            int(np.sum(res["valid"])),
            len(res["samples"]),
            int(np.sum(res["inside"])),
        ),
        file=out,
    )
    for name, pop in populations(res):
        summary[name] = dict(
            joints=int(np.sum(pop)), next_to_single=int(np.sum(pop & single))
        )
        for suffix, what in (("", "samples inside the masters"), ("_all", "all samples")):
            if suffix and np.all(res["inside"]):
                continue
            print(
                "-- %s: %d joints, %d of them next to a single quadratic; %s"
                % (name, int(np.sum(pop)), int(np.sum(pop & single)), what),
                file=out,
            )
            print(HEADER + "   (degrees)", file=out)
            rows = [
                ("cubic instances: kink", "kc", pop),
                ("quadratic instances: kink", "kq", pop),
                ("rounded quadratic: kink", "kr", pop),
                ("change by cu2qu (all)", "dq", pop),
                ("change by cu2qu (next to single)", "dq", pop & single),
                ("change by cu2qu (no single)", "dq", pop & ~single),
                ("  single + spline", "dq", pop & (one_single & other_spline)),
                ("  single + single", "dq", pop & both_single),
                ("  single + line", "dq", pop & (one_single & other_line)),
                ("change by rounding (all)", "dr", pop),
                ("change by rounding (next to single)", "dr", pop & single),
                ("no-single variant vs cubic", "dl", pop),
            ]
            sm = summary[name].setdefault("all samples" if suffix else "inside", {})
            for label, key, mask in rows:
                d = describe(res[key + suffix][mask])
                sm[label] = d
                print("%-36s %s" % (label, fmt(d)), file=out)
            print(HEADER + "   (font units, sideways)", file=out)
            for label, key, mask in [
                ("change by cu2qu (next to single)", "xq", pop & single),
                ("change by rounding (next to single)", "xr", pop & single),
                ("change by rounding (all)", "xr", pop),
            ]:
                d = describe(res[key + suffix][mask])
                sm[label + ", units"] = d
                print("%-36s %s" % (label, fmt(d)), file=out)
            m = pop & single
            if np.any(m):
                worse = int(np.sum(res["dq" + suffix][m] > res["dr" + suffix][m]))
                print(
                    "next to a single quadratic: cu2qu changes the angle more than rounding at %d of %d (%.2f%%)"
                    % (worse, int(np.sum(m)), 100.0 * worse / np.sum(m)),
                    file=out,
                )
                sm["cu2qu more than rounding"] = worse
            print(
                "largest difference between the no-single variant and the cubic instances: %.3g degrees; with no single quadratic at the joint: %.3g"
                % (
                    np.max(res["dl" + suffix][pop]) if np.any(pop) else 0,
                    np.max(res["dq" + suffix][pop & ~single]) if np.any(pop & ~single) else 0,
                ),
                file=out,
            )
    return summary


RULE_KS = (0.0, 0.01, 0.03, 0.1, 0.25, 0.5, 1.0)
GAP_KS = (0.1, 0.25, 0.5, 1.0, 2.0)


def rules(data, meta, base, out=sys.stdout):
    """Candidate rules: what they cost and what they leave."""
    kind = data["kind"]
    n = data["n"]
    prev = data["prev"]
    tol = meta["tol"]
    cubic = kind == 1
    single = cubic & (n == 1)
    points = int(len(kind) + np.sum(n[cubic]))
    dev = single_deviation(data)
    # Only the joints next to a single quadratic can change.
    only = single.copy()
    ok = prev >= 0
    only[ok] |= single[prev[ok]]
    sub = analyze(data, meta, None, only)
    gap = single_deviation(data, pair_fraction_gap)
    cands = [("as is", None)]
    for k in RULE_KS:
        cands.append(("deviation <= %g tol" % k, single & (dev > k * tol)))
    for k in GAP_KS:
        cands.append(("fractions <= %g tol" % k, single & (gap > k * tol)))
    cands.append(("no single quadratic", single.copy()))
    print(
        "-- candidate rules for single quadratics (tolerance %g, %d points, %d single of %d cubics)"
        % (tol, points, int(np.sum(single)), int(np.sum(cubic))),
        file=out,
    )
    summary = {}
    pop_names = [name for name, _ in populations(sub)]
    print(
        "%-22s %9s %8s   %s"
        % (
            "rule",
            "+points",
            "+%",
            "then, change by cu2qu at joints "
            + " | ".join("%s: >0.1 >0.5 >1 max(deg) max(units)" % nm for nm in pop_names),
        ),
        file=out,
    )
    for label, bump in cands:
        res = sub if bump is None else analyze(data, meta, bump, only)
        extra = 0 if bump is None else int(np.sum(bump))
        cells = []
        summary[label] = dict(extra_points=extra, percent=100.0 * extra / points)
        for name, pop in populations(res):
            x = res["dq"][pop]
            cell = (
                int(np.sum(x > 0.1)),
                int(np.sum(x > 0.5)),
                int(np.sum(x > 1.0)),
                float(np.max(x)) if len(x) else 0.0,
                float(np.max(res["xq"][pop])) if len(x) else 0.0,
            )
            summary[label][name] = cell
            cells.append("%7d %6d %6d %8.3f %7.3f" % cell)
        print(
            "%-22s %9d %7.3f%%   %s"
            % (label, extra, 100.0 * extra / points, " | ".join(cells)),
            file=out,
        )
    # A single quadratic that is linear in the control points instead.
    mid = analyze(data, meta, None, only, midpoint=True)
    print(
        "-- a linear single quadratic, off-curve point (3 (c1 + c2) - c0 - c3) / 4, for the same curves: change of the angle at joints, in the masters themselves and in instances",
        file=out,
    )
    summary["linear single"] = {}
    for name, pop in populations(mid):
        row = {}
        for label, res, key in (
            ("intersection, masters", sub, "dqm"),
            ("intersection, instances", sub, "dq"),
            ("linear, masters", mid, "dqm"),
            ("linear, instances", mid, "dq"),
        ):
            x = res[key][pop]
            cell = (
                int(np.sum(x > 0.1)),
                int(np.sum(x > 0.5)),
                int(np.sum(x > 1.0)),
                float(np.max(x)) if len(x) else 0.0,
            )
            row[label] = cell
            print(
                "%-32s %-26s joints %8d  >0.1 %7d  >0.5 %7d  >1 %7d  max %8.3f"
                % ((name, label, len(x)) + cell),
                file=out,
            )
        summary["linear single"][name] = row
    d = dev[single]
    if len(d):
        qs = np.percentile(d, (50, 90, 99, 99.9, 100))
        print(
            "deviation of the single quadratics, in units: p50 %.4f p90 %.4f p99 %.4f p99.9 %.4f max %.4f; exactly zero for %d of %d"
            % (tuple(qs) + (int(np.sum(d == 0)), len(d))),
            file=out,
        )
        summary["deviation"] = dict(
            zip(("p50", "p90", "p99", "p99.9", "max"), map(float, qs)),
            zero=int(np.sum(d == 0)),
            count=len(d),
        )
    return summary


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    path = args[0]
    data, meta, names = load(path)
    res = analyze(data, meta)
    n = data["n"]
    kind = data["kind"]
    cubics = int(np.sum(kind == 1))
    hist = {int(k): int(v) for k, v in zip(*np.unique(n[kind == 1], return_counts=True))}
    print(
        "font %s: %d masters, tolerance %g, %d glyphs, %d segments, %d cubics"
        % (meta["name"], data["pts"].shape[1], meta["tol"], meta["glyphs"], len(kind), cubics)
    )
    print("quadratics per cubic: %s" % hist)
    points = int(len(kind) + np.sum(n[kind == 1]))
    print("points (one on-curve per segment, n off-curve per cubic): %d" % points)
    summary = report(res, meta["name"])
    out = dict(
        meta=meta,
        masters=int(data["pts"].shape[1]),
        cubics=cubics,
        n_hist=hist,
        points=points,
        samples=len(res["samples"]),
        samples_inside=int(np.sum(res["inside"])),
        summary=summary,
    )
    if "--rules" in sys.argv:
        del res
        out["rules"] = rules(data, meta, None)
    if len(args) > 1:
        with open(args[1], "w") as f:
            json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
