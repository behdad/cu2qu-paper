#!/bin/bash
# The real master groups of kparity.sh, converted with fontTools upstream main
# (the base of the pull request) and with the branch: which groups differ?
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad/agent-kinks
A=/home/behdad/.cache/cu2qu-work/agent-kinks
D=$A/kparity-out
set -e
FT_TREE=$A/fonttools-main bash $S/py.sh $S/kparity_py.py run $D/inputs.txt $D/py-main.txt
wc -l $D/py-main.txt $D/py.txt
paste -d'\n' $D/py-main.txt $D/py.txt | awk '
  NR % 2 == 1 {a = $0; next}
  {
    total++
    if (a != $0) {
      differ++
      na = split(a, x, " "); nb = split($0, y, " ")
      if (na != nb) count++; else points++
    }
  }
  END {printf "groups %d, differ %d: %d in the number of points, %d in coordinates only\n", total, differ + 0, count + 0, points + 0}'
