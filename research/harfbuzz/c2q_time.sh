#!/bin/bash
# Usage: c2q_time.sh <worktree> <prototype-script|v1|none> <rounds>
# Builds a variant of src/hb-gpu-cu2qu.hh WITHOUT counters and runs the
# interleaved hb-gpu-encode-all A/B against the base libs ($BASE).
set -e
S=$(cd "$(dirname "$0")" && pwd)
WT=$1
FIX=$2
ROUNDS=${3:-5}
cd "$WT"
git checkout -q src/hb-gpu-cu2qu.hh
case "$FIX" in
  none) ;;
  v1) python3 "$S/c2q_fix_v1.py" src/hb-gpu-cu2qu.hh ;;
  *) python3 "$S/$FIX" src/hb-gpu-cu2qu.hh nostats ;;
esac
ninja -C build src/libharfbuzz-gpu.so.0.61520.0 2>&1 | grep -E " error|warning: " | head -8 || true
sed -i 's|^BASE=\$S/base/src-export/build/src$|BASE=${BASE:-$S/base/src-export/build/src}|' "$S/bench.sh"
bash "$S/bench.sh" "$WT" "$ROUNDS" \
  /usr/share/fonts/google-noto-sans-sc-fonts/NotoSansSC-Regular.otf=2 \
  /usr/share/fonts/adobe-source-code-pro-fonts/SourceCodePro-Regular.otf=100 \
  /usr/share/fonts/stix-fonts/STIXTwoMath-Regular.otf=15 \
  /home/behdad/.fonts/NotoSansCJKtc-VF.otf=1 \
  /usr/share/fonts/google-noto/NotoSans-Regular.ttf=50 2>&1 | tail -5
