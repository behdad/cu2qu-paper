#!/bin/bash
# Differential dumps of cu2qu results on the fixed input set (odd inputs left
# out, see ft_odd.py): compiled old, compiled new, pure old, pure new.
W=/home/behdad/hb/cu2qu-work
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
PY="$W/venv/bin/python"
export SKIP_ODD=1
run () { (cd /tmp && PYTHONPATH="$2" "$PY" "$S/ft_equiv.py" "$W/out.$1.txt" "$S/ft_cubics.json" 2>&1 | tail -2 | tr '\n' ' '; echo); }
run old-cy "$W/ft-old-cy/Lib" &
run new-cy "$W/ft-new-cy/Lib" &
run old-pure "$S/ft-old/Lib" &
run new-pure /home/behdad/fonttools/Lib &
wait
wc -l "$W"/out.*.txt | sed "s|$W/||"
cmpf () {
  if cmp -s "$W/out.$1.txt" "$W/out.$2.txt"; then echo "$1 == $2: IDENTICAL"
  else echo "$1 vs $2: differ in $(diff "$W/out.$1.txt" "$W/out.$2.txt" | grep -c '^<') lines"; fi
}
cmpf old-cy old-pure
cmpf new-cy new-pure
cmpf old-cy new-cy
cmpf old-pure new-pure
