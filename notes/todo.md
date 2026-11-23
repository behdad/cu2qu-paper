# To do

## Before the note goes on the site

- Behdad to read it. The title, the slug, the byline and the "About this
  note" paragraph are his to word.
- Pushing and deploying the homepage checkout are Behdad's. After a
  deployment, `curl -sI https://behdad.org/doc/cu2qu` must say `text/html`
  (not `text/javascript`), and the page must be looked at once in a browser.
- The entry in `doc/index.html` is a plain one, first under 2026 (the index
  is newest first; not featured, by Behdad's word), with the note's subtitle
  as its tagline: "The error of cu2qu's spline in closed form."
- The stylesheet link in the page is `../style.css?v=589bd39f2d`. The value
  is the first ten digits of `git hash-object style.css` in the homepage
  checkout; it changes when the stylesheet does. The figure script is linked
  the same way; `build.sh` sets that one, and `check_page.py` checks it.
- Section 8 says fontTools pull request 4452 and kurbo pull request 625 are
  "proposed". Change the wording when they are merged, and the two reference
  entries with them.

## Done after the note went live (2026-11-23)

- The note was deployed at https://behdad.org/doc/cu2qu. The first deployment
  had the script as `cu2qu.js`, which made the extensionless URL ambiguous
  (see `decisions.md`); renamed to `cu2qu-figures.js`, to be deployed again.
- At Behdad's word ("after we publish, put links in fonttools / kurbo /
  harfbuzz (merged) PRs"), one line was added to the description of each of
  fontTools 4452, kurbo 625 and HarfBuzz 6879, 6880, 6881 and 6882, pointing
  to the note (6880, 6881 and 6882 at the sections they use).

## Done after the two reviews (2026-11-23)

- See `reviews.md`. Two versions of the page followed the reviews; each
  needs deploying by Behdad.
- Not done: the figures have not been tried on a phone, only in a Chrome
  that emulates one (`scripts/test_touch.js`).

## Found while writing

- **fontTools and kurbo: a sentence in the documentation of the error factor
  was wrong for three to eight quadratics.** Acted on, 2026-11-23, at
  Behdad's "fix the comments". The docstring of `spline_error_factor` (and
  the doc comment of the kurbo port) said the polynomial of the pieces in
  between is "largest for k = 1, at v = sqrt(c/3)". That is so from n = 9.
  For fewer the largest value of the second piece is where it meets the
  first piece (v = -1/2), which is larger than the value at sqrt(c/3) and
  below the first piece's largest. The number returned is right for every
  n; only the sentence was off. Fixed by fontTools commit 41b05dbf2 in pull
  request 4452 and kurbo commit e44347f in pull request 625, comment only.
  The measurements here were made at 4201aba56 and 5f6ba66, one commit
  before each; the code is the same.
- Compiled with fused multiply-add (`gcc -mfma`), fontTools' cu2qu differs
  from pure Python in the last bit on non-integer input; see
  `notes/fonttools.md`. Not in the note.
- The page has no social card of its own; it uses `doc/og-card.png`.

## Could be added to the note

- A figure of the fit test that does not return (the piece that is its own
  half).
- The cost of forbidding single quadratics for masters (2.4% more points in
  the first kink measurement; not measured again).
- Rounding to integers: how often the rounded spline is over the tolerance
  (848 of 77,533 cubics at a tolerance of one unit in a research run, largest
  1.45; not measured again).
