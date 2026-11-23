#!/bin/bash
# Usage: run_c2q.sh <lib-dir> <font>...
# cu2qu counters (needs the temporary instrumentation) and blob sizes from
# util/hb-gpu-encode-all, per font.
S=$(cd "$(dirname "$0")" && pwd)
LIB=$1; shift
for font in "$@"; do
  echo "== $(basename "$font")"
  HB_CU2QU_STATS=1 LD_LIBRARY_PATH=$LIB "$S/encode_all" "$font" > "$S/c2q.out" 2> "$S/c2q.err"
  grep -E "^cu2qu|^  " "$S/c2q.err"
  grep -E "total:|avg:" "$S/c2q.out" | tr '\n' ' '; echo
done
