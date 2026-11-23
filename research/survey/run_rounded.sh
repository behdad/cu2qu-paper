#!/bin/bash
# cu2qu's splines with the off-curve points rounded to integers, as a glyf table
# stores them (implied on-curve points re-derived from the rounded off-curves).
set -e
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=$S/research-survey
B=/home/behdad/hb/cu2qu-work/research-survey
BIN=$B/target/release/c2q-survey
for tol in 1.0 0.5; do
  python3 $W/round_offcurves.py $B/res/fonttools_cu2qu.$tol.bin $B/rounded_fonttools_cu2qu.$tol.bin
  $BIN eval fonttools_cu2qu_rounded $tol $S/ft_cubics.json $B/rounded_fonttools_cu2qu.$tol.bin > $W/out/extra.fonttools_rounded.$tol.json
done
python3 - <<EOF
import json
for tol in ("1.0", "0.5"):
    d = json.load(open("$W/out/extra.fonttools_rounded.%s.json" % tol))
    print(
        "cu2qu at", tol, "then off-curves rounded:",
        "Hausdorff/tol mean %.3f p99 %.3f max %.3f;" % (d["haus"]["mean"], d["haus"]["p99"], d["haus"]["max"]),
        "cubics over tol", d["viol"], "(over 1.1 tol: %d);" % d["viol_10pct"],
        "end tangent off > 1 deg:", d["tan_gt_1"],
        "; tangent deg mean %.3f p99 %.2f" % (d["tan_deg"]["mean"], d["tan_deg"]["p99"]),
    )
EOF
