#!/bin/bash
# Usage: c2q_v3.sh <worktree> <prototype-script> [fonts...]
# Like c2q_variant.sh, for a prototype script that carries its own counters.
set -e
S=$(cd "$(dirname "$0")" && pwd)
WT=$1; shift
FIX=$1; shift
cd "$WT"
git checkout -q src/hb-gpu-cu2qu.hh
python3 "$S/$FIX" src/hb-gpu-cu2qu.hh
ninja -C build src/libharfbuzz-gpu.so.0.61520.0 2>&1 | grep -E " error|warning: " | head -8 || true
if [ $# -eq 0 ]; then
  set -- /usr/share/fonts/google-noto-sans-sc-fonts/NotoSansSC-Regular.otf \
         /usr/share/fonts/adobe-source-code-pro-fonts/SourceCodePro-Regular.otf \
         /usr/share/fonts/stix-fonts/STIXTwoMath-Regular.otf \
         /home/behdad/.fonts/NotoSansCJKtc-VF.otf
fi
bash "$S/run_c2q.sh" "$WT/build/src" "$@"
