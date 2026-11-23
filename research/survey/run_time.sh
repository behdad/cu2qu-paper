#!/bin/bash
# Single-thread conversion time per cubic (best of 7 passes over the corpus), Rust harness.
set -e
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=$S/research-survey
B=/home/behdad/hb/cu2qu-work/research-survey
BIN=$B/target/release/c2q-survey
for tol in 1.0 0.5 0.1; do
  for m in kurbo_to_quads lyon lyon017 kurbo_approx_spline reimpl_hb_direct reimpl_truong_presub reimpl_colomitchi_half; do
    $BIN time $m $tol $S/ft_cubics.json 7
  done
done | tee $W/out/time.rust.jsonl
