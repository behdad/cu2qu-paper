#!/bin/bash
# Do kurbo (the agent's clone, as checked out) and fontTools (the agent's
# clone, pure Python) convert real masters together to the same points, bit
# for bit?  usage: kparity.sh [LIMIT_PER_FONT]
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad/agent-kinks
A=/home/behdad/.cache/cu2qu-work/agent-kinks
D=$A/kparity-out
LIMIT=${1:-200000}
mkdir -p $D
set -e
bash $S/py.sh $S/kparity_py.py gen $D/inputs.txt $LIMIT $A/data/*.npz $A/data-tt/*.npz
(cd $A/kparity && CARGO_TARGET_DIR=$A/kparity-target cargo +1.99.0 build --release -q -j 8)
$A/kparity-target/release/kparity < $D/inputs.txt > $D/kurbo.txt &
bash $S/py.sh $S/kparity_py.py run $D/inputs.txt $D/py.txt &
wait
wc -l $D/inputs.txt $D/kurbo.txt $D/py.txt
if cmp -s $D/kurbo.txt $D/py.txt; then
  echo "kurbo and fontTools: identical output"
else
  echo "kurbo and fontTools differ on $(paste -d'\n' $D/kurbo.txt $D/py.txt | awk 'NR % 2 == 1 {a = $0; next} a != $0 {n++} END {print n + 0}') groups"
fi
