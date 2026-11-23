# To do

## Before the note goes on the site

- Behdad to read it. The title, the slug and the "About this note" paragraph
  are his to word.
- Put this repository on GitHub (Behdad). The page's Appendix B links to
  `https://github.com/behdad/cu2qu-paper` and to files under
  `blob/main/scripts/`; change the links if the name differs. Until it is
  there, `check_page.py --links --repository-pending` lets those links fail.
- The entry in `doc/index.html` is a plain one under 2026. The other papers
  of the year are `featured-work` entries with a tagline; make it one if it
  should be, for example with: "The error of cu2qu's spline in closed form:
  the number of quadratics is computed, not searched for."
- The stylesheet link in the page is `../style.css?v=589bd39f2d`. The value
  is the first ten digits of `git hash-object style.css` in the homepage
  checkout; it changes when the stylesheet does.
- Section 8 says fontTools pull request 4452 and kurbo pull request 625 are
  "proposed". Change the wording when they are merged, and the two reference
  entries with them.

## Found while writing, not acted on

- **fontTools and kurbo: a sentence in the documentation of the error factor
  is wrong for four to eight quadratics.** The docstring of
  `spline_error_factor` (and the doc comment of the kurbo port) says the
  polynomial of the pieces in between is "largest for k = 1, at
  v = sqrt(c/3)". That is so from n = 9. For n from 4 to 8 the largest value
  of the second piece is where it meets the first piece (v = -1/2), which is
  larger than the value at sqrt(c/3). The number returned is right for every
  n, because the first piece's maximum is larger than both up to n = 54; only
  the sentence is off. The note has it right (section 3). A comment-only fix
  for both pull requests.
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
