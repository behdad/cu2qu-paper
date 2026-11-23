#!/bin/bash
# Final check of the committed kurbo branch: bit-for-bit against the fontTools
# branch on both input sets, then the local CI checks.
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=/home/behdad/hb/cu2qu-work
D="$W/kdiff-out"
(cd "$S/kdiff/new" && CARGO_TARGET_DIR="$W/kdiff-target-new" cargo +1.99.0 build --release -q -j 8) || exit 1
"$W/kdiff-target-new/release/kdiff-new" < "$D/inputs.txt" > "$D/kurbo-new.txt" || echo "CRASH on the main set"
"$W/kdiff-target-new/release/kdiff-new" < "$D/single_ties.txt" > "$D/single_ties.kurbo-new.txt" || echo "CRASH on the single-quadratic limits"
cmp "$D/kurbo-new.txt" "$D/py-new.txt" && echo "main set: kurbo branch == fontTools branch, bit for bit ($(wc -l < "$D/py-new.txt") inputs)"
cmp "$D/single_ties.kurbo-new.txt" "$D/single_ties.py-new.txt" && echo "single-quadratic limits: kurbo branch == fontTools branch, bit for bit ($(wc -l < "$D/single_ties.py-new.txt") inputs)"
bash "$S/kurbo_checks.sh" 2>&1 | grep -E "^===|exit" | paste - - | sed 's/=== //' | awk -F'\t' '{printf "%-58s %s\n", $1, $2}'
