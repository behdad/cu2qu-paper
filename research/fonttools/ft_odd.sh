#!/bin/bash
# Outcomes on odd inputs, one process per input: compiled old against
# compiled new, and pure old against pure new.
W=/home/behdad/hb/cu2qu-work
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
PY="$W/venv/bin/python"
run () { (cd /tmp && PYTHONPATH="$2" "$PY" "$S/ft_odd.py" "$PY" > "$W/odd.$1.txt" 2>&1); }
run old-cy "$W/ft-old-cy/Lib" &
run new-cy "$W/ft-new-cy/Lib" &
run old-pure "$S/ft-old/Lib" &
run new-pure /home/behdad/fonttools/Lib &
wait
for v in old-cy new-cy old-pure new-pure; do echo "$v: $(tail -1 "$W/odd.$v.txt")"; done
echo "== compiled: old -> new, where they differ (outcome pairs, counted)"
paste -d'|' "$W/odd.old-cy.txt" "$W/odd.new-cy.txt" | awk -F'|' '$1 != $2' | grep -v "^# summary" | sed -E 's/^.{26} tol .{13} all_quadratic=(True |False) (curves|curve )  //' | awk -F'|' '{ sub(/^.{26} tol .{13} all_quadratic=(True |False) (curves|curve )  /, "", $2); print $1 "  ->  " $2 }' | sort | uniq -c
echo "== compiled new: inputs that still crash"
grep CRASH "$W/odd.new-cy.txt" | cut -c1-120 | head -20
echo "== pure vs compiled, new: lines that differ"
diff "$W/odd.new-pure.txt" "$W/odd.new-cy.txt" | grep -c '^<'
