#!/bin/sh
# Every run of the measuring instrument that the tables read.  About 20 minutes on 24 threads.
#   gcc -O2 -ffp-contract=off -fopenmp -o $W/gt gt.c -lm
R=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad/research-error
W=/home/behdad/hb/cu2qu-work/research-error
set -e
$W/gt $W/corpus.bin $W/gt_ramp.bin 16 0 0.1 64 256
$W/gt $W/corpus.bin $W/gt_zero.bin 16 1 0.1 64 256
$W/gt $W/corpus.bin $W/gt_alt.bin 16 2 0.1 64 256
sh $R/run_sets.sh 8 0 ramp
for n in A B M; do
  $W/gt $W/interp/$n.bin $W/interp/$n.gt 10 0 0 64 256
done
echo done
