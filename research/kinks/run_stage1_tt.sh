#!/bin/bash
# Stage 1 on cubic masters made from TrueType variable fonts (their quadratics
# raised to cubics, rounded to integers or not), into a separate data
# directory, and stage 2 on them.
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad/agent-kinks
D=/home/behdad/.cache/cu2qu-work/agent-kinks/data-tt
R=$S/results-tt
mkdir -p $D $R
one() {
  name=$1; kind=$2; path=$3
  bash $S/py.sh $S/stage1.py $name $kind "$path" -1 $D/$name.npz 8 > $D/$name.log 2>&1 || { echo "FAILED stage1 $name"; tail -3 $D/$name.log; return; }
  bash $S/py.sh $S/stage2.py $D/$name.npz $R/$name.json --rules > $R/$name.txt 2>&1 || echo "FAILED stage2 $name"
}
one tt-notosans ttvf "/usr/share/fonts/google-noto-vf/NotoSans[wght].ttf" &
one tt-notosans-exact ttvf-exact "/usr/share/fonts/google-noto-vf/NotoSans[wght].ttf" &
one tt-notosansitalic ttvf "/usr/share/fonts/google-noto-vf/NotoSans-Italic[wght].ttf" &
one tt-adwaitasans ttvf /usr/share/fonts/adwaita-sans-fonts/AdwaitaSans-Regular.ttf &
one tt-adwaitasansitalic ttvf /usr/share/fonts/adwaita-sans-fonts/AdwaitaSans-Italic.ttf &
one tt-bitter ttvf /home/behdad/.fonts/Bitter-VariableFont_wght.ttf &
one tt-vazirmatn ttvf "/usr/share/fonts/vazirmatn-vf-fonts/Vazirmatn[wght].ttf" &
one tt-notonaskharabic ttvf "/usr/share/fonts/google-noto-vf/NotoNaskhArabic[wght].ttf" &
wait
echo done
