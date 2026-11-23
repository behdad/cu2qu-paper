#!/bin/sh
# Check the note and put its figures in, in the homepage checkout:
# doc/cu2qu.html and doc/cu2qu-figures.js.
#
# Usage: sh build.sh [SITE]        (SITE defaults to ~/behdad.org)
#
# It checks the formulas, compares the figure script with them, writes the
# figures' static copies into the page, and checks the page's tables and
# quoted numbers against data/. If there is a google-chrome, it then lays
# the page out at phone and print widths, works the figures' controls, and
# drags their control points by touch. It stops at the first check that
# fails. It does not rerun the measurements (see README.md) and commits
# nothing.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
SITE=${1:-$HOME/behdad.org}
PAGE=$SITE/doc/cu2qu.html
JS=$SITE/doc/cu2qu-figures.js

cd "$HERE/scripts"
python3 check_math.py > "$HERE/data/check_math.txt"
python3 check_js.py "$JS" > "$HERE/data/check_js.txt"
node render_figures.js "$PAGE" "$JS" "$HERE/data/figures.json"
# the page links the figure script by the first ten digits of its git hash,
# as the site does for its stylesheet, so that a changed script is fetched
VERSION=$(git -C "$SITE" hash-object "$JS" | cut -c1-10)
sed -i "s|cu2qu-figures.js?v=[0-9a-f]*\"|cu2qu-figures.js?v=$VERSION\"|" "$PAGE"
python3 check_page.py "$PAGE" "$HERE/data"
if command -v google-chrome > /dev/null; then
  python3 check_layout.py "$PAGE"
  python3 test_interaction.py "$PAGE" > "$HERE/data/interaction.txt"
  node test_touch.js "$SITE"
fi
