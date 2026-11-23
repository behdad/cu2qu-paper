#!/bin/bash
# Run the new cu2qu tests against the OLD fontTools (upstream main export), to
# see which of them the old code fails. The old code has no
# spline_error_factor, so the test file gets a copy of it.
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
D="$S/old-tests"
rm -rf "$D"; mkdir -p "$D/Tests/cu2qu"
cp -r /home/behdad/fonttools/Tests/cu2qu/data "$D/Tests/cu2qu/"
python3 - "$D/Tests/cu2qu/cu2qu_test.py" <<'EOF'
import sys
s = open('/home/behdad/fonttools/Tests/cu2qu/cu2qu_test.py').read()
a = "from fontTools.cu2qu.cu2qu import MAX_N, spline_error_factor\n"
assert a in s
s = s.replace(a, '''from fontTools.cu2qu.cu2qu import MAX_N


def spline_error_factor(n):
    g = (n - 2) / (4 * (n - 1))
    error = (4 / 27) * (1 + g) * (1 + g) * (1 + g)
    if n > 3:
        c = (n - 3) / (4 * (n - 1))
        w = c / 3
        error = max(error, c + 2 * w * math.sqrt(w))
    return error
''')
open(sys.argv[1], 'w').write(s)
EOF
cd "$D" && PYTHONPATH="${OLD_LIB:-$S/ft-old/Lib}" "${PYTHON:-python3}" -m pytest Tests/cu2qu/cu2qu_test.py -q -p no:cacheprovider -p no:randomly 2>&1 | grep -E "^(FAILED|ERROR)|passed|failed|error" | cut -c1-200
