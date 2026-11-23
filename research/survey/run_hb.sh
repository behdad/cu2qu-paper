#!/bin/bash
# HarfBuzz hb_gpu_cubic_to_quadratics() (unmodified header copy + shim) at the three tolerances.
set -e
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=$S/research-survey
B=/home/behdad/hb/cu2qu-work/research-survey
for tol in 1.0 0.5 0.1; do
  $B/bin/hb_c2q $B/corpus.bin $tol $B/res/harfbuzz.$tol.bin 5
done | tee $W/out/time.harfbuzz.jsonl
