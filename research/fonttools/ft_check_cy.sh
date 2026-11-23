#!/bin/bash
# Compiled-module checks for the cu2qu change:
#  1. Cython's own warnings for cu2qu.py, old against new
#  2. Tests/cu2qu on the compiled new tree, and on the pure-Python working tree
#  3. the new tests against the compiled OLD tree
#  4. differential dumps: compiled old, compiled new, pure new
W=/home/behdad/hb/cu2qu-work
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
PY="$W/venv/bin/python"
echo "== Cython warnings for cu2qu.py (old / new)"
grep -c "^warning: .*cu2qu/cu2qu.py" "$W/ft-old-cy.build.log"
grep -c "^warning: .*cu2qu/cu2qu.py" "$W/ft-new-cy.build.log"
grep "^warning: .*cu2qu/cu2qu.py" "$W/ft-new-cy.build.log" | head -5
echo "== Tests/cu2qu, compiled new"
(cd "$W/ft-new-cy" && PYTHONPATH="$W/ft-new-cy/Lib" "$PY" -c "from fontTools.cu2qu.cu2qu import COMPILED; assert COMPILED" && PYTHONPATH="$W/ft-new-cy/Lib" "$PY" -m pytest Tests/cu2qu -q -p no:cacheprovider 2>&1 | tail -2)
echo "== Tests/cu2qu, pure new (working tree)"
(cd /home/behdad/fonttools && PYTHONPATH=/home/behdad/fonttools/Lib "$PY" -c "from fontTools.cu2qu.cu2qu import COMPILED; assert not COMPILED" && PYTHONPATH=/home/behdad/fonttools/Lib "$PY" -m pytest Tests/cu2qu -q -p no:cacheprovider 2>&1 | tail -2)
echo "== new tests against compiled OLD (expect the tie and NaN tests to fail; NaN may crash)"
OLD_LIB="$W/ft-old-cy/Lib" PYTHON="$PY" bash "$S/ft_old_tests.sh" 2>&1 | tail -6
echo "exit status of that run: $?"
if [ "$1" = dumps ]; then
  echo "== differential dumps"
  (cd /tmp && PYTHONPATH="$W/ft-old-cy/Lib" "$PY" "$S/ft_equiv.py" "$W/out.old-cy.txt" "$S/ft_cubics.json" 2>&1 | tail -2) &
  (cd /tmp && PYTHONPATH="$W/ft-new-cy/Lib" "$PY" "$S/ft_equiv.py" "$W/out.new-cy.txt" "$S/ft_cubics.json" 2>&1 | tail -2) &
  (cd /tmp && PYTHONPATH=/home/behdad/fonttools/Lib "$PY" "$S/ft_equiv.py" "$W/out.new-pure.txt" "$S/ft_cubics.json" 2>&1 | tail -2) &
  wait
  wc -l "$W"/out.*.txt
  cmp -s "$W/out.new-cy.txt" "$W/out.new-pure.txt" && echo "new: compiled == pure, IDENTICAL" || echo "new: compiled vs pure DIFFER in $(diff "$W/out.new-cy.txt" "$W/out.new-pure.txt" | grep -c '^<') lines"
  cmp -s "$W/out.old-cy.txt" "$W/out.new-cy.txt" && echo "compiled: old == new, IDENTICAL" || echo "compiled: old vs new differ in $(diff "$W/out.old-cy.txt" "$W/out.new-cy.txt" | grep -c '^<') lines"
fi
