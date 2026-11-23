#!/bin/bash
# The whole fontTools test suite (as tox runs it: Tests + doctests), on the
# working tree in pure Python and on the compiled export of it.
# Usage: ft_fulltests.sh [pure|cy|old-cy|all]
W=/home/behdad/hb/cu2qu-work
PY="$W/venv/bin/python"
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$W/tmp"
mkdir -p "$TMPDIR"
run () {
  name=$1; tree=$2; want=$3
  (
    cd "$tree" || exit 1
    PYTHONPATH="$tree/Lib" "$PY" -c "from fontTools.cu2qu.cu2qu import COMPILED; import fontTools, sys; assert COMPILED == $want, COMPILED; print('$name: fontTools from', fontTools.__file__, 'compiled', COMPILED)"
    rm -rf "$W/pytest-tmp-$name"
    PYTHONPATH="$tree/Lib" "$PY" -m pytest -q -p no:cacheprovider --basetemp="$W/pytest-tmp-$name" > "$W/fulltests.$name.log" 2>&1
    echo "$name: exit $?: $(tail -1 "$W/fulltests.$name.log")"
    grep -E "^(FAILED|ERROR) " "$W/fulltests.$name.log" | head -20
  )
}
which=${1:-all}
if [ "$which" = pure ] || [ "$which" = all ]; then run pure /home/behdad/fonttools False & fi
if [ "$which" = cy ] || [ "$which" = all ]; then run cy "$W/ft-new-cy" True & fi
if [ "$which" = old-cy ]; then run old-cy "$W/ft-old-cy" True & fi
wait
