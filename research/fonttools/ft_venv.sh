#!/bin/bash
# Make a virtualenv with Cython, the test tools and fontTools' test
# requirements (those that install from wheels), on the home disk: the
# scratchpad is on a tmpfs with a quota.
W=/home/behdad/hb/cu2qu-work
set -e
mkdir -p "$W"
python3 -m venv "$W/venv"
"$W/venv/bin/pip" install -q --upgrade pip
"$W/venv/bin/pip" install -q cython setuptools wheel pytest pytest-randomly black==26.3.1 mypy
set +e
# One at a time, so that one without a wheel for this Python does not stop the rest.
grep -v '^#' /home/behdad/fonttools/requirements.txt | sed 's/ *#.*//' | grep -v '^$' | while read -r req; do
  if "$W/venv/bin/pip" install -q --only-binary=:all: "$req" >/dev/null 2>&1; then
    echo "ok      $req"
  else
    echo "FAILED  $req"
  fi
done
"$W/venv/bin/python" -c "import Cython, pytest, black; print('cython', Cython.__version__, 'pytest', pytest.__version__, 'black', black.__version__)"
"$W/venv/bin/pip" list 2>/dev/null | grep -i -E "^(fonttools|ufoLib2|ufo2ft|glyphsLib|lxml|scipy|sympy|skia-pathops|uharfbuzz|freetype-py|brotli|zopfli|fs|unicodedata2) "
du -sh "$W/venv"
