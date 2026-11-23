#!/bin/bash
# A/B of the single-quadratic skip in hb-gpu's cu2qu: build main's header,
# keep its libs as the base, build the changed header, then run the encoder
# benchmark (instruction counts, interleaved).
WT=/home/behdad/harfbuzz/.claude/worktrees/gpu-cu2qu-timeout
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=/home/behdad/.cache/cu2qu-work
set -e
cd "$WT"
cp src/hb-gpu-cu2qu.hh "$S/hb-gpu-cu2qu.new.hh"
git show origin/main:src/hb-gpu-cu2qu.hh > src/hb-gpu-cu2qu.hh
ninja -C build > "$S/hb_build_base.log" 2>&1 || { tail -20 "$S/hb_build_base.log"; cp "$S/hb-gpu-cu2qu.new.hh" src/hb-gpu-cu2qu.hh; exit 1; }
rm -rf "$W/hb-base"; mkdir -p "$W/hb-base"
cp -a build/src/libharfbuzz*.so* "$W/hb-base/"
cp "$S/hb-gpu-cu2qu.new.hh" src/hb-gpu-cu2qu.hh
ninja -C build > "$S/hb_build_new.log" 2>&1 || { tail -20 "$S/hb_build_new.log"; exit 1; }
grep -c "warning" "$S/hb_build_new.log" | sed 's/^/warnings in the build of the change: /'
git status --short | grep -v '^??'
mkdir -p "$S/hbbench"
cp /home/behdad/.claude/projects/-home-behdad-harfbuzz/memory/gpu-tools/bench.sh "$S/hbbench/bench.sh"
BASE="$W/hb-base" bash "$S/hbbench/bench.sh" "$WT" 5 \
  /usr/share/fonts/google-noto-sans-sc-fonts/NotoSansSC-Regular.otf=2 \
  /usr/share/fonts/adobe-source-code-pro-fonts/SourceCodePro-Regular.otf=40 \
  /usr/share/fonts/stix-fonts/STIXTwoMath-Regular.otf=15 \
  /usr/share/fonts/google-noto-sans-cjk-vf-fonts/NotoSansCJK-VF.ttc=1 \
  ${CONTROL_FONT:-/usr/share/fonts/google-noto-vf/NotoSans[wght].ttf}=10
