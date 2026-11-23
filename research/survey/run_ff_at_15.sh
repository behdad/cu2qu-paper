#!/bin/bash
# FontForge's own tolerance is 1.5 units when it rounds control points to integers,
# so compare it with fontTools cu2qu run at max_err = 1.5.
set -e
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=$S/research-survey
B=/home/behdad/hb/cu2qu-work/research-survey
BIN=$B/target/release/c2q-survey
PYTHONPATH=$S/ft-old/Lib python3 $W/harness/ft_c2q.py $S/ft_cubics.json 1.5 $B/res15_fonttools_cu2qu.bin > /dev/null
$BIN gen kurbo_to_quads 1.5 $S/ft_cubics.json $B/res15_kurbo_to_quads.bin
$BIN eval fonttools_cu2qu_at_1.5 1.5 $S/ft_cubics.json $B/res15_fonttools_cu2qu.bin > $W/out/extra.fonttools_at_1.5.json
$BIN eval kurbo_to_quads_at_1.5 1.5 $S/ft_cubics.json $B/res15_kurbo_to_quads.bin > $W/out/extra.kurbo_to_quads_at_1.5.json
$BIN eval fontforge_vs_1.5 1.5 $S/ft_cubics.json $B/res/fontforge.1.0.bin > $W/out/extra.fontforge_vs_1.5.json
python3 - <<EOF
import json
for f in ("extra.fonttools_at_1.5", "extra.kurbo_to_quads_at_1.5", "extra.fontforge_vs_1.5"):
    d = json.load(open("$W/out/" + f + ".json"))
    print(
        d["label"],
        "quads/cubic %.4f" % d["mean_n"],
        "TT/cubic %.4f" % d["tt_mean"],
        "H/tol mean %.3f p99 %.3f max %.3f" % (d["haus"]["mean"], d["haus"]["p99"], d["haus"]["max"]),
        "over tol", d["viol"],
        "tan>1deg", d["tan_gt_1"],
    )
EOF
