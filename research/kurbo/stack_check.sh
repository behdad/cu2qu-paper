#!/bin/bash
# The stacked kurbo branch (~/kurbo working tree: halving + computed n) against
# the fontTools branch (~/fonttools working tree), every coordinate by its
# bits, on the main input set and the single-quadratic limits; then the local
# CI checks and the instruction count.
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=/home/behdad/.cache/cu2qu-work
D="$W/kdiff-out"
if [ ! -s "$D/py-new.txt" ]; then
  (cd /tmp && PYTHONPATH=/home/behdad/fonttools/Lib python3 "$S/kdiff_py.py" run "$D/inputs.txt" "$D/py-new.txt" 2>&1 | tail -1) &
  (cd /tmp && PYTHONPATH=/home/behdad/fonttools/Lib python3 "$S/kdiff_py.py" run "$D/single_ties.txt" "$D/single_ties.py-new.txt" 2>&1 | tail -1) &
  wait
fi
(cd "$S/kdiff/new" && CARGO_TARGET_DIR="$W/kdiff-target-new" cargo +1.99.0 build --release -q -j 8) || exit 1
"$W/kdiff-target-new/release/kdiff-new" < "$D/inputs.txt" > "$D/kurbo-stack.txt" || echo "CRASH on the main set"
"$W/kdiff-target-new/release/kdiff-new" < "$D/single_ties.txt" > "$D/single_ties.kurbo-stack.txt" || echo "CRASH on the single-quadratic limits"
echo "== main set, kurbo branch vs fontTools branch"
python3 "$S/kdiff/analyze.py" "$D/inputs.txt" "$D/kurbo-stack.txt" "$D/py-new.txt" | cut -c1-250
echo "== single-quadratic limits (accuracy within 3 ulps of the single quadratic's error)"
paste -d'\n' "$D/single_ties.kurbo-stack.txt" "$D/single_ties.py-new.txt" | awk 'NR % 2 == 1 {a = $0; next} {t++; if (a != $0) n++} END {print "lines", t, "differ", n + 0}'
echo "== local CI checks (failures only)"
bash "$S/kurbo_checks.sh" 2>&1 | grep -E "^===|exit" | paste - - | sed 's/=== //' | awk -F'\t' '{printf "%-58s %s\n", $1, $2}' | grep -v "exit 0"
echo "== instructions"
bash "$S/kdiff/bench.sh" stack /home/behdad/kurbo/kurbo 2>&1 | tail -1
