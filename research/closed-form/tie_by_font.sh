#!/bin/bash
# Usage: tie_by_font.sh <instrumented-lib-dir> <font-list>
# Which fonts produce the "trial n != formula n" lines, how many, and how
# close to the tolerance they are.
S=$(cd "$(dirname "$0")" && pwd)
LIB=$1
LIST=$2
total=0
while read -r font; do
  out=$(LD_LIBRARY_PATH=$LIB MAX_GLYPHS=400 "$S/gpu_dump" "$font" 2>&1 >/dev/null | grep "^trial n")
  [ -z "$out" ] && continue
  n=$(echo "$out" | wc -l)
  total=$((total + n))
  ratios=$(echo "$out" | awk '{print $12}' | sort -u | head -4 | tr '\n' ' ')
  kinds=$(echo "$out" | awk '{print $3 $6}' | sort | uniq -c | tr '\n' ' ')
  echo "$n $(basename "$font") | $kinds| ratios: $ratios"
done < "$LIST" | sort -rn | head -14
