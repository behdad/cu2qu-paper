#!/bin/bash
# Export two fontTools trees to the home-disk work dir and build their Cython
# extensions in place:
#   ft-old-cy: upstream main (origin/main)
#   ft-new-cy: the cu2qu-compute-n working tree (HEAD plus uncommitted changes
#              to the files this branch touches)
# Usage: ft_build_cy.sh [old|new|both]
W=/home/behdad/hb/cu2qu-work
SRC=/home/behdad/fonttools
which=${1:-both}
set -e
build () {
  name=$1; rev=$2
  rm -rf "$W/$name"; mkdir -p "$W/$name"
  git -C "$SRC" archive "$rev" | tar -x -C "$W/$name"
  if [ "$name" = ft-new-cy ]; then
    cp "$SRC/Lib/fontTools/cu2qu/cu2qu.py" "$W/$name/Lib/fontTools/cu2qu/cu2qu.py"
    cp "$SRC/Tests/cu2qu/cu2qu_test.py" "$W/$name/Tests/cu2qu/cu2qu_test.py"
  fi
  (cd "$W/$name" && FONTTOOLS_WITH_CYTHON=1 "$W/venv/bin/python" setup.py -q build_ext --inplace > "$W/$name.build.log" 2>&1) || { tail -30 "$W/$name.build.log"; exit 1; }
  grep -c -i "warning" "$W/$name.build.log" | sed "s/^/$name: warnings in build log: /"
  grep -i "cu2qu" "$W/$name.build.log" | grep -i "warning" | head -5
  ls "$W/$name/Lib/fontTools/cu2qu/" | grep -E "\.so$"
  (cd /tmp && PYTHONPATH="$W/$name/Lib" "$W/venv/bin/python" -c "
import fontTools.cu2qu.cu2qu as m
print('$name', m.__file__.replace('$W/', ''), 'COMPILED', m.COMPILED)
assert m.COMPILED")
}
if [ "$which" = old ] || [ "$which" = both ]; then build ft-old-cy origin/main; fi
if [ "$which" = new ] || [ "$which" = both ]; then build ft-new-cy HEAD; fi
