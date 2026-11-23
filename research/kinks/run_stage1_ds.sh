#!/bin/bash
# Stage 1 on designspaces and .glyphs sources of noto-source, in parallel.
# usage: run_stage1_ds.sh NAME:PATH ...
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad/agent-kinks
D=/home/behdad/.cache/cu2qu-work/agent-kinks/data
mkdir -p $D
for spec in "$@"; do
  name=${spec%%:*}
  path=${spec#*:}
  (bash $S/py.sh $S/stage1.py $name ds $path -1 $D/$name.npz 1 > $D/$name.log 2>&1; echo "$name: $(tail -c 400 $D/$name.log)") &
done
wait
