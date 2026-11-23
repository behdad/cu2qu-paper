#!/bin/bash
# Benchmark cu2qu: old against new, compiled and pure, interleaved, three
# rounds each (each number is already a best-of-7 inside ft_bench.py).
W=/home/behdad/hb/cu2qu-work
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
PY="$W/venv/bin/python"
for round in 1 2 3; do
  for v in old-cy new-cy old-pure new-pure; do
    case $v in
      old-cy) lib="$W/ft-old-cy/Lib";;
      new-cy) lib="$W/ft-new-cy/Lib";;
      old-pure) lib="$S/ft-old/Lib";;
      new-pure) lib=/home/behdad/fonttools/Lib;;
    esac
    echo "$v $(cd /tmp && PYTHONPATH="$lib" "$PY" "$S/ft_bench.py" "$S/ft_cubics.json")"
  done
done
echo "== fontTools' own benchmark (python -m fontTools.cu2qu.benchmark)"
for v in old-cy new-cy; do
  case $v in
    old-cy) lib="$W/ft-old-cy/Lib";;
    new-cy) lib="$W/ft-new-cy/Lib";;
  esac
  echo "$v: $(cd /tmp && PYTHONPATH="$lib" "$PY" -m fontTools.cu2qu.benchmark | tr '\n' ' ')"
done
