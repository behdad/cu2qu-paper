#!/bin/sh
# Run the measuring instrument on every generated set, ramp weighting, n = 1..8, no early stop.
W=/home/behdad/hb/cu2qu-work/research-error
NMAX=${1:-8}
MODE=${2:-0}
TAG=${3:-ramp}
DIR=${4:-sets}
for f in $W/$DIR/*.bin; do
  name=$(basename $f .bin)
  $W/gt $f $W/$DIR/$name.gt_$TAG $NMAX $MODE 0 64 512
done
