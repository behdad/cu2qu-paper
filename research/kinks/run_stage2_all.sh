#!/bin/bash
# Stage 2 (with candidate rules) on every data set of stage 1.
# usage: run_stage2_all.sh [DATA_DIR [RESULTS_DIR]]
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad/agent-kinks
D=${1:-/home/behdad/.cache/cu2qu-work/agent-kinks/data}
R=${2:-$S/results}
mkdir -p $R
run() {
  name=$(basename $1 .npz)
  bash $S/py.sh $S/stage2.py $1 $R/$name.json --rules > $R/$name.txt 2>&1 || echo "FAILED $name"
}
# The two CJK fonts take several GB each: one at a time, the rest in parallel.
for f in $D/*.npz; do
  case $f in
    *cjk*) ;;
    *) run $f & ;;
  esac
done
wait
for f in $D/*cjk*.npz; do
  [ -e "$f" ] && run $f
done
echo done
