#!/bin/sh
# Check the note and put its figures in, in the homepage checkout:
# doc/cu2qu.html and doc/cu2qu.js.
#
# Usage: sh build.sh [SITE]        (SITE defaults to ~/behdad.org)
#
# It checks the formulas, compares the figure script with them, writes the
# figures' static copies into the page, and checks the page's tables and
# quoted numbers against data/. It stops at the first check that fails. It
# does not rerun the measurements (see README.md) and commits nothing.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
SITE=${1:-$HOME/behdad.org}
PAGE=$SITE/doc/cu2qu.html
JS=$SITE/doc/cu2qu.js

cd "$HERE/scripts"
python3 check_math.py > "$HERE/data/check_math.txt"
python3 check_js.py "$JS" > "$HERE/data/check_js.txt"
node render_figures.js "$PAGE" "$JS" "$HERE/data/figures.json"
python3 check_page.py "$PAGE" "$HERE/data"
