#!/bin/bash
# Verify the depth limit in fontTools' fit test (working tree of ~/fonttools)
# against the branch's HEAD commit:
#   black; Tests/cu2qu pure and compiled; results on the fixed 515k input set
#   identical to HEAD (compiled and pure); the endless-recursion inputs, one
#   process each, on compiled HEAD and compiled new; whole suite pure and
#   compiled.
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=/home/behdad/.cache/cu2qu-work
SRC=/home/behdad/fonttools
PY="$W/venv/bin/python"
ulimit -c 0
echo "=== black"
(cd "$SRC" && "$PY" -m black --check Lib/fontTools/cu2qu/cu2qu.py Tests/cu2qu/cu2qu_test.py 2>&1 | tail -1)
echo "=== build HEAD and working tree with Cython"
rm -rf "$W/ft-head-cy"; mkdir -p "$W/ft-head-cy"
git -C "$SRC" archive HEAD | tar -x -C "$W/ft-head-cy"
(cd "$W/ft-head-cy" && FONTTOOLS_WITH_CYTHON=1 "$PY" setup.py -q build_ext --inplace > "$W/ft-head-cy.build.log" 2>&1) || { tail -5 "$W/ft-head-cy.build.log"; exit 1; }
bash "$S/ft_build_cy.sh" new 2>&1 | tail -1
grep -c "^warning: .*cu2qu/cu2qu.py" "$W/ft-new-cy.build.log" | sed 's/^/Cython warnings for cu2qu.py: /'
echo "=== Tests/cu2qu"
(cd "$W/ft-new-cy" && PYTHONPATH="$W/ft-new-cy/Lib" "$PY" -m pytest Tests/cu2qu fontTools.cu2qu -q -p no:cacheprovider 2>&1 | tail -1 | sed 's/^/compiled: /')
(cd "$SRC" && PYTHONPATH="$SRC/Lib" "$PY" -m pytest Tests/cu2qu fontTools.cu2qu -q -p no:cacheprovider 2>&1 | tail -1 | sed 's/^/pure:     /')
echo "=== the inputs that recursed without end, one process each"
for v in head-cy new-cy; do
  for call in "[(0,0),(0,0),(0,0),(25,0)], 100/27" "[(-45,-3),(-30,-2),(0,5),(0,9)], 20/3" "[(0,0),(0,0),(0,0),(1282,0)], 189.92592592592595"; do
    out=$(cd /tmp && PYTHONPATH="$W/ft-$v/Lib" "$PY" -c "
from fontTools.cu2qu import cu2qu
from fontTools.cu2qu.cu2qu import curve_to_quadratic
assert cu2qu.COMPILED
print(len(curve_to_quadratic($call)) - 2, 'quadratic(s)')" 2>/dev/null)
    echo "$v  curve_to_quadratic($call): ${out:-CRASH (exit $?)}"
  done
done
echo "=== results on the fixed input set: HEAD compiled, new compiled, new pure"
export SKIP_ODD=1
run () { (cd /tmp && PYTHONPATH="$2" "$PY" "$S/ft_equiv.py" "$W/out.$1.txt" "$S/ft_cubics.json" > /dev/null 2>&1); }
run head-cy "$W/ft-head-cy/Lib" &
run new-cy "$W/ft-new-cy/Lib" &
run new-pure "$SRC/Lib" &
wait
wc -l < "$W/out.new-cy.txt" | sed 's/^/inputs: /'
cmp -s "$W/out.head-cy.txt" "$W/out.new-cy.txt" && echo "HEAD == new (compiled): IDENTICAL" || echo "HEAD vs new (compiled): DIFFER"
cmp -s "$W/out.new-cy.txt" "$W/out.new-pure.txt" && echo "new compiled == new pure: IDENTICAL" || echo "new compiled vs pure: DIFFER"
rm -f "$W"/out.head-cy.txt "$W"/out.new-cy.txt "$W"/out.new-pure.txt
echo "=== whole suite"
bash "$S/ft_fulltests.sh" all 2>&1 | tail -2
