#!/bin/bash
# Everything, in order.  Needs: the built Rust harness (cargo build --release with
# CARGO_TARGET_DIR=$B/target), build_cpp.sh, node, python3 with fontforge.
#   run_all.sh [quick]     quick = skip the external generators (results already there)
set -e
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=$S/research-survey
B=/home/behdad/hb/cu2qu-work/research-survey
mkdir -p $B/res $B/percubic $W/out

if [ "$1" != "quick" ]; then
  bash $W/run_rust.sh                 # kurbo, lyon, reimplementations
  bash $W/run_ft.sh                   # fontTools (pure Python)
  bash $W/run_js.sh                   # cubic2quad (node)
  bash $W/run_hb.sh                   # HarfBuzz header
  bash $W/run_skia.sh                 # Skia GrPathUtils
  python3 $W/harness/ff_c2q.py $S/ft_cubics.json $B/res/fontforge.1.0.bin > $W/out/time.fontforge.jsonl
fi
METHODS="oracle_ramp oracle_qm oracle_truong oracle_tangent oracle_tangent_implied" bash $W/run_rust.sh
bash $W/run_eval.sh
bash $W/run_time.sh > /dev/null
python3 $W/summarize.py > $W/out/summary.txt
python3 $W/analyze_percubic.py > $W/out/analysis.txt
$B/target/release/c2q-survey consts > $W/out/consts.txt
$B/target/release/c2q-survey samen $S/ft_cubics.json > $W/out/samen.txt
$B/target/release/c2q-survey rate > $W/out/rate.txt
echo all done
