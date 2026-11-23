#!/bin/bash
# fontTools cu2qu (pure Python from the clean upstream export) at the three tolerances.
set -e
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=$S/research-survey
B=/home/behdad/hb/cu2qu-work/research-survey
for tol in 1.0 0.5 0.1; do
  PYTHONPATH=$S/ft-old/Lib python3 $W/harness/ft_c2q.py $S/ft_cubics.json $tol $B/res/fonttools_cu2qu.$tol.bin
done | tee $W/out/time.fonttools.jsonl
