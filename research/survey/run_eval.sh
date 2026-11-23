#!/bin/bash
# Evaluate every results file with the common evaluator.
#   run_eval.sh [method ...]      (default: every *.bin in $B/res)
set -e
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=$S/research-survey
B=/home/behdad/hb/cu2qu-work/research-survey
BIN=$B/target/release/c2q-survey
CORPUS=$S/ft_cubics.json
mkdir -p $W/out $B/percubic

if [ $# -gt 0 ]; then
  FILES=""
  for m in "$@"; do
    for tol in 1.0 0.5 0.1; do
      [ -f $B/res/$m.$tol.bin ] && FILES="$FILES $B/res/$m.$tol.bin"
    done
  done
else
  FILES=$(ls $B/res/*.bin)
fi

for f in $FILES; do
  base=$(basename $f .bin)       # method.tol
  tol=${base##*.}
  tol=$(echo $base | sed -E 's/.*\.([0-9]+\.[0-9]+)$/\1/')
  m=$(echo $base | sed -E 's/\.([0-9]+\.[0-9]+)$//')
  $BIN eval $m $tol $CORPUS $f $B/percubic/$base.txt > $W/out/eval.$base.json
  echo "$base ok"
done
