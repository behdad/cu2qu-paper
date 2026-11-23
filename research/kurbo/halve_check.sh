#!/bin/bash
# The standalone halving branch of kurbo (~/kurbo working tree) against
# fontTools MAIN, bit for bit, on the main input set and on the
# single-quadratic limits; plus upstream kurbo against the same for contrast.
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=/home/behdad/.cache/cu2qu-work
D="$W/kdiff-out"
mkdir -p "$D"
[ -s "$D/inputs.txt" ] || python3 "$S/kdiff_py.py" gen "$D/inputs.txt" "$S/ft_cubics.json" 2>&1 | tail -1
[ -s "$D/single_ties.txt" ] || python3 "$S/kdiff/gen_single_ties.py" "$D/single_ties.txt" 40000 2>&1 | tail -1
if [ ! -s "$D/py-old.txt" ]; then
  (cd /tmp && PYTHONPATH="$S/ft-old/Lib" python3 "$S/kdiff_py.py" run "$D/inputs.txt" "$D/py-old.txt" 2>&1 | tail -1) &
  (cd /tmp && PYTHONPATH="$S/ft-old/Lib" python3 "$S/kdiff_py.py" run "$D/single_ties.txt" "$D/single_ties.py-old.txt" 2>&1 | tail -1) &
  wait
fi
(cd "$S/kdiff/new" && CARGO_TARGET_DIR="$W/kdiff-target-new" cargo +1.99.0 build --release -q -j 8) || exit 1
"$W/kdiff-target-new/release/kdiff-new" < "$D/inputs.txt" > "$D/kurbo-halve.txt" || echo "CRASH on the main set"
"$W/kdiff-target-new/release/kdiff-new" < "$D/single_ties.txt" > "$D/single_ties.kurbo-halve.txt" || echo "CRASH on the single-quadratic limits"
wc -l "$D/py-old.txt" "$D/kurbo-halve.txt" "$D/single_ties.py-old.txt" "$D/single_ties.kurbo-halve.txt" | sed "s|$D/||"
cmp "$D/kurbo-halve.txt" "$D/py-old.txt" && echo "main set: kurbo (halving fix) == fontTools main, bit for bit"
cmp "$D/single_ties.kurbo-halve.txt" "$D/single_ties.py-old.txt" && echo "single-quadratic limits: kurbo (halving fix) == fontTools main, bit for bit"
python3 "$S/kdiff/analyze.py" "$D/inputs.txt" "$D/kurbo-halve.txt" "$D/py-old.txt" | cut -c1-260
