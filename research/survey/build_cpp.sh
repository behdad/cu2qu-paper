#!/bin/bash
# Build the two C++ harnesses.
#  - hb_c2q:   an unmodified copy of HarfBuzz's src/hb-gpu-cu2qu.hh (copied here,
#              the worktree is not touched) compiled against harness/hb/ (shim).
#  - skia_c2q: links Skia's own GrPathUtils from the prebuilt libskia.a.
set -e
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=$S/research-survey
B=/home/behdad/hb/cu2qu-work/research-survey
HB=/home/behdad/harfbuzz/.claude/worktrees/gpu-cu2qu-timeout/src
mkdir -p $B/bin $B/hbcopy

cp $HB/hb-gpu-cu2qu.hh $B/hbcopy/hb-gpu-cu2qu.hh
cmp $HB/hb-gpu-cu2qu.hh $B/hbcopy/hb-gpu-cu2qu.hh
c++ -std=c++17 -O2 -I$W/harness/hb -I$B/hbcopy -o $B/bin/hb_c2q $W/harness/hb_c2q.cpp

c++ -std=c++20 -O2 -fno-exceptions -fno-rtti -DNDEBUG -DSK_GANESH -DSK_GL -DSK_R32_SHIFT=16 \
    -DSK_GAMMA_APPLY_TO_A8 -DSK_ALLOW_STATIC_GLOBAL_INITIALIZERS=1 -DGPU_TEST_UTILS=1 \
    -I/home/behdad/skia -o $B/bin/skia_c2q $W/harness/skia_c2q.cpp \
    /home/behdad/skia/out/Release/libskia.a -lpthread -ldl
echo built
