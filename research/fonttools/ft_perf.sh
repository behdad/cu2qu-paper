#!/bin/bash
# Count user-space instructions of a fixed cu2qu workload for several builds:
#   old-cy, head-cy (the branch's HEAD commit), new-cy (working tree),
#   old-pure, head-pure, new-pure.
# Prints instructions per call (workload minus the load-only baseline).
W=/home/behdad/hb/cu2qu-work
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
PY="$W/venv/bin/python"
count () {  # lib mode -> instructions
  (cd /tmp && PYTHONPATH="$1" perf stat -x, -e instructions:u "$PY" "$S/ft_perf.py" "$S/ft_cubics.json" "$2" 2>&1 >/dev/null | awk -F, '/instructions/ {print $1}')
}
calls () { case $1 in single|half) echo 232599;; multi) echo 77532;; random) echo 60000;; esac; }
printf "%-10s %12s %12s %12s %12s\n" build single half multi random
for v in ${VARIANTS:-old-cy head-cy new-cy old-pure head-pure new-pure}; do
  case $v in
    old-cy) lib="$W/ft-old-cy/Lib";;
    head-cy) lib="$W/ft-head-cy/Lib";;
    new-cy) lib="$W/ft-new-cy/Lib";;
    old-pure) lib="$S/ft-old/Lib";;
    head-pure) lib="$W/ft-head-pure/Lib";;
    new-pure) lib=/home/behdad/fonttools/Lib;;
  esac
  base=$(count "$lib" none)
  line="$v"
  for mode in single half multi random; do
    n=$(count "$lib" $mode)
    line="$line $(( (n - base) / $(calls $mode) ))"
  done
  printf "%-10s %12s %12s %12s %12s\n" $line
done
