#!/bin/bash
# cubic2quad (JavaScript, the downloaded index.js run unmodified) at the three tolerances.
set -e
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=$S/research-survey
B=/home/behdad/hb/cu2qu-work/research-survey
for tol in 1.0 0.5 0.1; do
  node $W/harness/c2q_js.js $B/dl/cubic2quad/index.js $S/ft_cubics.json $tol $B/res/cubic2quad_js.$tol.bin
done | tee $W/out/time.cubic2quad.jsonl
