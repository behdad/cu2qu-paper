#!/bin/bash
# Independent brute-force check of the evaluator's Hausdorff distance on 2000
# pseudo-randomly chosen cubics per results file (dense sampling of both curves,
# 2000 points per segment, point-to-polyline distances, no solver, no refinement).
# (Not run on reimpl_alfke_penner: its worst outputs have ~1000 segments and the
# brute force is quadratic in that.)
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=$S/research-survey
B=/home/behdad/hb/cu2qu-work/research-survey
BIN=$B/target/release/c2q-survey
: > $W/out/verify.txt
for f in fonttools_cu2qu.1.0 fonttools_cu2qu.0.1 harfbuzz.0.5 kurbo_to_quads.1.0 kurbo_to_quads.0.1 lyon.0.5 lyon017.0.5 skia.1.0 skia_cal.0.1 cubic2quad_js.1.0 cubic2quad_js.0.1 fontforge.1.0 reimpl_truong_presub.0.5 oracle_ramp.1.0 oracle_tangent.0.5; do
  $BIN verify $S/ft_cubics.json $B/res/$f.bin 2000 | sed "s|$B/res/||" >> $W/out/verify.txt
done
cat $W/out/verify.txt
