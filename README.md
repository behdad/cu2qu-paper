# How many quadratics does a cubic take?

What is behind the note at https://behdad.org/doc/cu2qu: the scripts that
check its mathematics and measure its numbers, their output, the tools that
build and check the page, the notes kept while working on it, and the
earlier research tools.

The note itself is not here. It lives in the homepage's repository, as
`doc/cu2qu.html` with its figure script `doc/cu2qu-figures.js`.

## Layout

| | |
|---|---|
| `scripts/` | The formulas in Python, the checks, the measurements, and the tools for the page. Listed below. |
| `data/` | The scripts' output as quoted in the note, and `cubics.json`, the 77,533 cubics of three fonts that the measurements run on. |
| `notes/` | What was decided and why, the research behind the note, what is left to do. Working notes: what the note itself says was checked again; much of the rest was not. |
| `research/` | Earlier tools, from the work on HarfBuzz, fontTools and kurbo that led here and from the research runs. Kept as they were: paths in them point to where they ran. |
| `build.sh` | Runs the checks and writes the figures' static copies into the page. |

## The scripts

The formulas

- `cu2qu_math.py`: the formulas of the note in plain Python: the spline, the
  error function, the error factor and the limits, the count, the bound for
  one quadratic.
- `check_math.py`: checks every mathematical statement of the note against
  them, on random cubics and random weights, against fine grids, and in exact
  rational arithmetic where the statement is an identity of polynomials.
  Output: `data/check_math.txt`.

The fonts

- `corpus.py`: writes the cubics of some CFF fonts to a JSON file.
  `data/cubics.fonts.txt` names the three fonts and their versions.
- `measure.py`: the tables and numbers about those cubics, from the formulas
  alone; with `--hausdorff` also the distance between each cubic and its
  splines as sets of points. Output: `data/measure.txt`.

fontTools, before and after

- `compare_revs.py`: runs two revisions of fontTools' cu2qu on the cubics and
  compares their results with each other and with the formula.
  Output: `data/compare-fonttools.txt`.
- `perf_fonttools.py`: instructions and time per call for the two revisions,
  compiled with Cython and in pure Python; uses `cu2qu_workload.py` and
  `perf`. Output: `data/perf-fonttools.txt`.
- `check_endless.py`: the inputs on which the fit test did not return.
  Output: `data/endless.txt`.

The revisions compared are f3f799e1b (main) and 4201aba56
(https://github.com/fonttools/fonttools/pull/4452).

Masters

- `kinks.py`: smooth joints of interpolated instances, before and after
  conversion, on the masters of 24 variable-font sources. Its input is those
  masters as arrays, which are not here (272 MB); see `notes/kinks.md`.
  Output: `data/kinks.txt`.

The page

- `render_figures.js`: runs the page's figure script under node and writes
  the figures, as they are before anything is moved, into the page, for
  readers without scripts. Writes `data/figures.json`.
- `check_js.py`: compares the figure script's formulas with `cu2qu_math.py`.
  Output: `data/check_js.txt`.
- `check_page.py`: checks that the page is well formed and that its tables
  and quoted numbers are the ones in `data/`; with `--links`, that every
  link answers.
- `test_interaction.py`: moves every slider and drags points in a headless
  Chrome, and prints what the figures report.

They need Python 3 with numpy and fontTools; node for the figures; Cython
and `perf` for `perf_fonttools.py`; google-chrome for `test_interaction.py`.

## Working on the note

Edit `doc/cu2qu.html` or `doc/cu2qu-figures.js` in the homepage checkout, then

    sh build.sh              # for ~/behdad.org; or: sh build.sh SITE

To look at it, serve the homepage checkout:

    python3 -m http.server 8000 --bind 127.0.0.1 --directory ~/behdad.org
    # http://localhost:8000/doc/cu2qu.html

Other checks, not part of `build.sh`:

    python3 scripts/check_page.py ~/behdad.org/doc/cu2qu.html data --links
    python3 scripts/test_interaction.py ~/behdad.org/doc/cu2qu.html

## Rerunning the measurements

From `scripts/`, with the paths of the machine the note was written on:

    python3 corpus.py ../data/cubics.json 3000 \
        /usr/share/fonts/adobe-source-code-pro-fonts/SourceCodePro-Regular.otf \
        /usr/share/fonts/stix-fonts/STIXTwoMath-Regular.otf \
        /usr/share/fonts/google-noto-sans-sc-fonts/NotoSansSC-Regular.otf
    python3 measure.py ../data/cubics.json --hausdorff > ../data/measure.txt
    python3 compare_revs.py ~/fonttools f3f799e1b 4201aba56 ../data/cubics.json \
        > ../data/compare-fonttools.txt
    python3 perf_fonttools.py ~/fonttools f3f799e1b 4201aba56 ../data/cubics.json \
        WORKDIR PYTHON-WITH-CYTHON > ../data/perf-fonttools.txt
    python3 check_endless.py ~/fonttools f3f799e1b 4201aba56 > ../data/endless.txt
    PYTHONPATH=~/fonttools/Lib python3 kinks.py MASTERS-DIR > ../data/kinks.txt

`measure.py --hausdorff` takes a few minutes; `perf_fonttools.py` builds
fontTools twice with Cython. If a number changes, `check_page.py` says which
sentence or table cell of the page no longer agrees.
