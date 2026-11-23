#!/bin/bash
# Usage: run_gdiff2.sh <worktree> <base-lib-dir> <tag>
# gpu_dump over the test fonts plus a few system CFF fonts, base vs new,
# compared with gdiff_c2q.py.
set -e
S=$(cd "$(dirname "$0")" && pwd)
WT=$1
BASE=$2
TAG=$3
FIX=$WT/build/src

g++ -O2 -o "$S/gpu_dump" "$S/gpu_dump.cc" -I"$WT/src" -L"$FIX" -lharfbuzz-gpu -lharfbuzz -lm

cd "$WT"
find test/api/fonts test/shape/data test/subset/data/fonts test/fuzzing/fonts -type f \
  \( -name '*.ttf' -o -name '*.otf' -o -name '*.ttc' -o -name '*.otc' -o -path 'test/fuzzing/fonts/*' \) \
  | sort > "$S/fonts.$TAG.txt"
for f in /usr/share/fonts/adobe-source-code-pro-fonts/SourceCodePro-Regular.otf \
	 /usr/share/fonts/stix-fonts/STIXTwoMath-Regular.otf \
	 /usr/share/fonts/google-noto-sans-sc-fonts/NotoSansSC-Regular.otf \
	 /usr/share/fonts/google-noto/NotoSans-Regular.ttf; do
  echo "$f" >> "$S/fonts.$TAG.txt"
done

export SYNTH=${SYNTH:-0}
export MAX_GLYPHS=${MAX_GLYPHS:-400}
( LD_LIBRARY_PATH=$BASE "$S/gpu_dump" $(cat "$S/fonts.$TAG.txt") > "$S/out.base.$TAG.txt" ) &
( LD_LIBRARY_PATH=$FIX  "$S/gpu_dump" $(cat "$S/fonts.$TAG.txt") > "$S/out.fix.$TAG.txt" ) &
wait
python3 "$S/gdiff_c2q.py" "$S" "$TAG"
