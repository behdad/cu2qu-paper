#!/bin/bash
# Convert the three input sets (general, ties, adversarial) with one
# fontTools tree, bit-exact results to $A/out/<label>.<set>.txt.
# Usage: run_ft.sh <label> <fontTools Lib dir> [set...]
W=/home/behdad/.cache/cu2qu-work
A=$W/agent-fit
SA=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad/agent-fit
label=$1; lib=$2; shift; shift
sets=${*:-general ties adversarial}
mkdir -p "$A/out"
cd /tmp || exit 1
for s in $sets; do
  ( PYTHONPATH="$lib" PYTHONDONTWRITEBYTECODE=1 "$W/venv/bin/python" "$SA/ft_run.py" "$A/data/inputs_$s.txt" "$A/out/$label.$s.txt" 2> "$A/out/$label.$s.log" ) &
done
wait
for s in $sets; do echo "$label $s: $(tail -1 "$A/out/$label.$s.log") $(head -1 "$A/out/$label.$s.log")"; done
