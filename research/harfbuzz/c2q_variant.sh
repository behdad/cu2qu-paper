#!/bin/bash
# Usage: c2q_variant.sh <worktree> <fix-script-or-none> [fonts...]
# Restores src/hb-gpu-cu2qu.hh, applies a candidate fix script (or none),
# then the counters, rebuilds the gpu library and prints the counters and
# blob sizes for the fonts.  Leaves the instrumented file in place.
set -e
S=$(cd "$(dirname "$0")" && pwd)
WT=$1; shift
FIX=$1; shift
cd "$WT"
git checkout -q src/hb-gpu-cu2qu.hh
if [ "$FIX" != none ]; then python3 "$S/$FIX" src/hb-gpu-cu2qu.hh; fi
python3 "$S/cu2qu_stats.py" src/hb-gpu-cu2qu.hh
ninja -C build src/libharfbuzz-gpu.so.0.61520.0 2>&1 | grep -E " error|warning: " | head -5 || true
if [ $# -eq 0 ]; then
  set -- /usr/share/fonts/google-noto-sans-sc-fonts/NotoSansSC-Regular.otf \
         /usr/share/fonts/adobe-source-code-pro-fonts/SourceCodePro-Regular.otf \
         /usr/share/fonts/stix-fonts/STIXTwoMath-Regular.otf \
         /home/behdad/.fonts/NotoSansCJKtc-VF.otf
fi
bash "$S/run_c2q.sh" "$WT/build/src" "$@"
