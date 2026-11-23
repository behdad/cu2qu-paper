#!/usr/bin/env python3
"""Tables from out/eval.<method>.<tol>.json (written by run_eval.sh)."""
import glob
import json
import os
import sys

W = os.path.dirname(os.path.abspath(__file__))
rows = {}
for f in glob.glob(os.path.join(W, "out", "eval.*.json")):
    d = json.load(open(f))
    rows[(d["label"], d["tol"])] = d

ORDER = [
    "fonttools_cu2qu",
    "kurbo_approx_spline",
    "kurbo_cubics_to_quadratic_splines",
    "harfbuzz",
    "reimpl_hb_direct",
    "reimpl_cu2qu_trial",
    "kurbo_to_quads",
    "lyon",
    "lyon017",
    "lyon017_monotonic",
    "skia",
    "skia_cal",
    "cubic2quad_js",
    "fontforge",
    "reimpl_colomitchi_half",
    "reimpl_colomitchi_adaptive",
    "reimpl_colomitchi_symmetric",
    "reimpl_truong_presub",
    "reimpl_alfke_penner",
    "oracle_ramp",
    "oracle_truong",
    "oracle_qm",
    "oracle_tangent_implied",
    "oracle_tangent",
]
labels = sorted({k[0] for k in rows}, key=lambda l: (ORDER.index(l) if l in ORDER else 99, l))
only = sys.argv[1:] or None

for tol in (1.0, 0.5, 0.1):
    print(f"\n=== tolerance {tol} font units; {next(iter(rows.values()))['cubics']} cubics ===")
    print(
        f"{'method':34s} {'quads/cu':>8s} {'TTpts/cu':>8s} {'expl.on':>7s} "
        f"{'H mean':>6s} {'H p50':>6s} {'H p99':>6s} {'H max':>7s} {'>tol':>6s} {'>1.1tol':>7s} "
        f"{'eq max':>7s} {'tan mean':>8s} {'tan p99':>7s} {'tan>1deg':>8s} {'kink max':>8s} {'kink>1':>6s} "
        f"{'bbox mean':>9s} {'bbox max':>8s}  n-histogram"
    )
    for l in labels:
        if only and l not in only:
            continue
        d = rows.get((l, tol))
        if not d:
            continue
        hist = d["hist"]
        hs = " ".join(f"{k}:{v}" for k, v in sorted(hist.items(), key=lambda kv: int(kv[0])) if int(k) <= 8)
        more = sum(v for k, v in hist.items() if int(k) > 8)
        if more:
            hs += f" >8:{more}"
        print(
            f"{l:34s} {d['mean_n']:8.4f} {d['tt_mean']:8.4f} {d['explicit_on']:7d} "
            f"{d['haus']['mean']:6.3f} {d['haus']['p50']:6.3f} {d['haus']['p99']:6.3f} {d['haus']['max']:7.3f} "
            f"{d['viol']:6d} {d['viol_10pct']:7d} "
            f"{d['eq']['max']:7.3f} {d['tan_deg']['mean']:8.3f} {d['tan_deg']['p99']:7.2f} {d['tan_gt_1']:8d} "
            f"{d['kink_deg']['max']:8.2f} {d['kink_gt_1']:6d} "
            f"{d['bbox']['mean']:9.4f} {d['bbox']['max']:8.3f}  {hs}"
        )
