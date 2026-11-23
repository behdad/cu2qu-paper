#!/bin/bash
# Bit-exact comparison of kurbo's cu2qu port with fontTools, old and new:
#   kurbo old  = upstream main (pristine clone)
#   kurbo new  = ~/kurbo working tree (branch cu2qu-compute-n)
#   py old     = fontTools upstream main
#   py new     = ~/fonttools working tree (branch cu2qu-compute-n)
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=/home/behdad/hb/cu2qu-work
D="$W/kdiff-out"
mkdir -p "$D"
set -e
if [ ! -s "$D/inputs.txt" ] || [ "$1" = regen ]; then
  python3 "$S/kdiff_py.py" gen "$D/inputs.txt" "$S/ft_cubics.json"
fi
(cd "$S/kdiff/old" && CARGO_TARGET_DIR="$W/kdiff-target-old" cargo +1.99.0 build --release -q -j 8)
(cd "$S/kdiff/new" && CARGO_TARGET_DIR="$W/kdiff-target-new" cargo +1.99.0 build --release -q -j 8)
"$W/kdiff-target-old/release/kdiff-old" < "$D/inputs.txt" > "$D/kurbo-old.txt" &
"$W/kdiff-target-new/release/kdiff-new" < "$D/inputs.txt" > "$D/kurbo-new.txt" &
(cd /tmp && PYTHONPATH="$S/ft-old/Lib" python3 "$S/kdiff_py.py" run "$D/inputs.txt" "$D/py-old.txt" 2>&1 | tail -1) &
(cd /tmp && PYTHONPATH=/home/behdad/fonttools/Lib python3 "$S/kdiff_py.py" run "$D/inputs.txt" "$D/py-new.txt" 2>&1 | tail -1) &
wait
wc -l "$D"/*.txt | sed "s|$D/||"
cmpf () {
  n=$(paste -d'\n' "$D/$1.txt" "$D/$2.txt" | awk 'NR % 2 == 1 {a = $0; next} a != $0 {n++} END {print n + 0}')
  echo "$1 vs $2: $n lines differ"
}
cmpf kurbo-old py-old
cmpf kurbo-new py-new
cmpf kurbo-old kurbo-new
cmpf py-old py-new
