#!/bin/bash
# Generate results with the Rust harness (real kurbo / lyon code and the
# labelled reimplementations) for the three tolerances.
set -e
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=$S/research-survey
B=/home/behdad/hb/cu2qu-work/research-survey
BIN=$B/target/release/c2q-survey
CORPUS=$S/ft_cubics.json

[ -f $B/corpus.bin ] || $BIN corpus2bin $CORPUS $B/corpus.bin

METHODS="${METHODS:-kurbo_to_quads kurbo_approx_spline kurbo_cubics_to_quadratic_splines lyon lyon017 lyon017_monotonic reimpl_hb_direct reimpl_cu2qu_trial reimpl_colomitchi_half reimpl_colomitchi_adaptive reimpl_colomitchi_symmetric reimpl_truong_presub reimpl_alfke_penner}"
TOLS="${TOLS:-1.0 0.5 0.1}"

for tol in $TOLS; do
  for m in $METHODS; do
    $BIN gen $m $tol $CORPUS $B/res/$m.$tol.bin
  done
done
echo done
