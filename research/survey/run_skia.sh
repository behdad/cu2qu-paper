#!/bin/bash
# Skia GrPathUtils::convertCubicToQuads, the real function from the prebuilt
# /home/behdad/skia/out/Release/libskia.a (see harness/skia_c2q.cpp and build_skia.sh).
#   skia            tolScale = tol          (the parameter taken at face value)
#   skia_cal        tolScale = 3.375 * tol  (calibrated: its acceptance test
#                   |c0 - c1| = |d3 of the piece| / 2 < tolScale bounds the
#                   equal-parameter error by (4/27) * 2 * tolScale = tol)
set -e
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=$S/research-survey
B=/home/behdad/hb/cu2qu-work/research-survey
for tol in 1.0 0.5 0.1; do
  $B/bin/skia_c2q $B/corpus.bin $tol $B/res/skia.$tol.bin 5
  cal=$(python3 -c "print($tol*3.375)")
  $B/bin/skia_c2q $B/corpus.bin $cal $B/res/skia_cal.$tol.bin 5
done | tee $W/out/time.skia.jsonl
