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
- `best_weights.py`: the best weights of section 4, by linear programming
  (needs scipy). Output: `data/best-weights.txt`.

The fonts

- `corpus.py`: writes the cubics of some CFF fonts to a JSON file.
  `data/cubics.fonts.txt` names the three fonts and their versions.
- `measure.py`: the tables and numbers about those cubics, from the formulas
  alone; with `--hausdorff` also the distance between each cubic and its
  splines as sets of points. Output: `data/measure.txt`.

fontTools, before and after

- `compare_revs.py`: runs two revisions of fontTools' cu2qu on the cubics and
  compares their results with each other and with the formulas: the number
  of quadratics, and the closed form of the single quadratic's error.
  Output: `data/compare-fonttools.txt`.
- `compare_revs_masters.py`: the same comparison on groups of compatible
  masters, converted together. Output: `data/compare-fonttools-masters.txt`.
- `ties.py`: what the recursive test did with splines exactly at a limit:
  accepted, rejected, or did not return. Output: `data/ties.txt`.
- `perf_fonttools.py`: instructions and time per call for the two revisions,
  compiled with Cython and in pure Python; uses `cu2qu_workload.py` and
  `perf`. Output: `data/perf-fonttools.txt`.
- `check_endless.py`: the inputs on which the fit test did not return.
  Output: `data/endless.txt`.

The revisions compared are f3f799e1b (main) and 4201aba56
(https://github.com/fonttools/fonttools/pull/4452).

kurbo, before and after

- `check_endless_kurbo.py`: the same inputs, and kurbo's own, in
  `CubicBez::approx_spline`. Output: `data/endless-kurbo.txt`.
- `compare_kurbo.py`: kurbo's results against fontTools', bit for bit.
  Output: `data/compare-kurbo.txt`.
- `perf_kurbo.py`: instructions per call. Output: `data/perf-kurbo.txt`.

The revisions are a5f97ff40 (main) and 5f6ba6679
(https://github.com/linebender/kurbo/pull/625). The three scripts export
them from a kurbo checkout and build small programs against them with cargo.

Masters

- `kinks.py`: smooth joints of interpolated instances, before and after
  conversion, on the masters of 24 variable fonts. Its input is those
  masters as arrays, which are not here (272 MB); see `notes/kinks.md`.
  Output: `data/kinks.txt`, and with `--smooth-within 1`,
  `data/kinks-within-one-degree.txt`.
- `kinks_inputs.py`: where those arrays came from, with versions.
  Output: `data/kinks-inputs.txt`.

The page

- `render_figures.js`: runs the page's figure script under node and writes
  the figures, as they are before anything is moved, into the page, for
  readers without scripts. Writes `data/figures.json`.
- `check_js.py`: compares the figure script's formulas with `cu2qu_math.py`.
  Output: `data/check_js.txt`.
- `check_page.py`: checks that the page is well formed and that its tables
  and quoted numbers are the ones in `data/`; with `--links`, that every
  link answers with a page.
- `check_layout.py`: lays the page out at phone, print and desktop widths in
  a headless Chrome, and looks for what does not fit.
- `test_interaction.py`: moves every slider, drags points and presses the
  reset buttons in a headless Chrome, with events it makes up, and checks
  what the figures report. Output: `data/interaction.txt`.
- `test_touch.js`: drags the control points with real touches, in a Chrome
  that emulates a phone, and with a real mouse.

They need Python 3 with numpy and fontTools; scipy for `best_weights.py`;
node for the figures (22 or later for `test_touch.js`); Cython and `perf`
for `perf_fonttools.py`; cargo and `perf` for the kurbo scripts;
google-chrome for the last three.

## Working on the note

Edit `doc/cu2qu.html` or `doc/cu2qu-figures.js` in the homepage checkout, then

    sh build.sh              # for ~/behdad.org; or: sh build.sh SITE

To look at it, serve the homepage checkout:

    python3 -m http.server 8000 --bind 127.0.0.1 --directory ~/behdad.org
    # http://localhost:8000/doc/cu2qu.html

`build.sh` also runs the three browser checks if there is a google-chrome.
Not part of it, because it asks other people's servers:

    python3 scripts/check_page.py ~/behdad.org/doc/cu2qu.html data --links

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
    python3 ties.py ~/fonttools f3f799e1b 4201aba56 > ../data/ties.txt
    python3 best_weights.py --table > ../data/best-weights.txt

    python3 check_endless_kurbo.py ~/kurbo a5f97ff 5f6ba66 WORKDIR > ../data/endless-kurbo.txt
    python3 compare_kurbo.py ~/kurbo 5f6ba66 ~/fonttools 4201aba56 ../data/cubics.json \
        WORKDIR > ../data/compare-kurbo.txt
    python3 perf_kurbo.py ~/kurbo a5f97ff 5f6ba66 ../data/cubics.json WORKDIR \
        > ../data/perf-kurbo.txt

    PYTHONPATH=~/fonttools/Lib python3 kinks.py MASTERS-DIR > ../data/kinks.txt
    PYTHONPATH=~/fonttools/Lib python3 kinks.py MASTERS-DIR --smooth-within 1 \
        > ../data/kinks-within-one-degree.txt
    python3 kinks_inputs.py MASTERS-DIR > ../data/kinks-inputs.txt
    python3 compare_revs_masters.py ~/fonttools f3f799e1b 4201aba56 MASTERS-DIR \
        > ../data/compare-fonttools-masters.txt

`measure.py --hausdorff` takes about twenty minutes; `perf_fonttools.py`
builds fontTools twice with Cython. If a number changes, `check_page.py` says
which sentence or table cell of the page no longer agrees.
