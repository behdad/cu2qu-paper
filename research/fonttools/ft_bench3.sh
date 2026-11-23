#!/bin/bash
# Benchmark three states of cu2qu, compiled and pure, interleaved:
#   old  = upstream main
#   head = the branch's HEAD commit (exported, without working-tree changes)
#   new  = the working tree
W=/home/behdad/hb/cu2qu-work
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
SRC=/home/behdad/fonttools
PY="$W/venv/bin/python"
set -e
# head: compiled and pure exports
rm -rf "$W/ft-head-cy" "$W/ft-head-pure"; mkdir -p "$W/ft-head-cy" "$W/ft-head-pure"
git -C "$SRC" archive HEAD | tar -x -C "$W/ft-head-cy"
git -C "$SRC" archive HEAD Lib/fontTools | tar -x -C "$W/ft-head-pure"
(cd "$W/ft-head-cy" && FONTTOOLS_WITH_CYTHON=1 "$PY" setup.py -q build_ext --inplace > "$W/ft-head-cy.build.log" 2>&1)
bash "$S/ft_build_cy.sh" new | tail -1
set +e
for round in 1 2 3; do
  for v in old-cy head-cy new-cy old-pure head-pure new-pure; do
    case $v in
      old-cy) lib="$W/ft-old-cy/Lib";;
      head-cy) lib="$W/ft-head-cy/Lib";;
      new-cy) lib="$W/ft-new-cy/Lib";;
      old-pure) lib="$S/ft-old/Lib";;
      head-pure) lib="$W/ft-head-pure/Lib";;
      new-pure) lib="$SRC/Lib";;
    esac
    echo "$v $(cd /tmp && PYTHONPATH="$lib" "$PY" "$S/ft_bench.py" "$S/ft_cubics.json")"
  done
done
echo "== fontTools' own benchmark (python -m fontTools.cu2qu.benchmark)"
for v in old-cy head-cy new-cy; do
  echo "$v: $(cd /tmp && PYTHONPATH="$W/ft-$v/Lib" "$PY" -m fontTools.cu2qu.benchmark | tr '\n' ' ')"
done
