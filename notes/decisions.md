# Decisions

What was decided while writing the note, by whom, and why. Dates are 2026.

## The page

- **Plain HTML with the site's stylesheet.** Behdad asked for simple HTML
  like the rest of the technical writing on his homepage, and then to just
  include his stylesheet. The page links `../style.css` and adds one small
  style block for what the stylesheet has no rules for (formulas, figures,
  tables, listings). It follows the hand-written pages of the site
  (`<body class="ltr">`, `site-nav`, `h1.title`, anchors in headings), not the
  Google Docs exports.
- **Relative links to the stylesheet and scripts.** The first draft linked
  `/style.css`, which only resolves from the site root and not when the file
  is opened directly. `../style.css` works from `/doc/cu2qu` and from a file.
- **MathJax.** The first draft had native MathML generated from TeX, set in
  Noto Sans Math. Behdad found it ugly and asked for MathJax. The formulas
  are now TeX in the page, typeset by MathJax 3.2.2 from jsDelivr. Its size is
  fixed (`matchFontHeight: false`, `scale: 1.06`): left to measure the text
  font at load, it came out 2.4 times too large in a headless Chrome, before
  the web fonts had arrived. Its menu is off, so that 300 formulas are not 300
  tab stops.
- **Live figures.** Behdad asked for interactive figures. One file,
  `cu2qu.js`, draws all five and wires sliders and dragging; no library. It
  is the only implementation: the copies in the page for readers without
  scripts are made by running the same file under node
  (`render_figures.js`), and `check_js.py` compares its formulas with the
  Python ones. An earlier Python figure script was dropped for this.
- **Figure colours.** The site's red for cu2qu's spline, and two more hues
  (blue `#2a78d6`, green `#1baf7a`) that pass a colour-blindness check against
  each other and the red on white. They fall back to fixed values when the
  stylesheet is missing, so that the figures never lose their curves.
- **Where things live.** Behdad, in the end: the results go in the homepage
  checkout and may be committed there; the research and reproduction notes
  and code go in this repository, which he puts on GitHub himself. So the
  page and `cu2qu.js` are only in the homepage checkout, this repository has
  no copy, and the page's Appendix B links here
  (`github.com/behdad/cu2qu-paper`). For a while the scripts were to be
  published next to the page instead (`doc/cu2qu-files/`); that was dropped.

## Who wrote it

- Claude wrote the note; Behdad directed the work. When told the index entry
  was a plain one that he might want to feature, he answered that the index
  is in reverse chronological order, and that the note is not to be featured
  "since you wrote it, not me". He gave it its tagline in the index: "The
  math behind cu2qu's algorithm". So: the entry is plain and first under 2026;
  the byline reads "Written by Claude, with Behdad Esfahbod"; "About this
  note" says the derivations, code, measurements, figures and text are
  Claude's, at his direction; and the few sentences that had him speaking in
  the first person are now in the third person or the passive. The byline
  was Claude's choice after that remark, not something he asked for in so
  many words; it is his to change.

## The content

- **Every number has a script.** Tables and quoted numbers are checked against
  `data/` by `check_page.py`. Numbers that come from the pull requests and
  were not measured again here are said to be the pull requests': the 773,396
  groups of masters, kurbo's 2,426 to 735 instructions, HarfBuzz's 14 to 19%.
- **Measurements were redone, not copied.** The corpus numbers of the research
  runs were reproduced by `measure.py`, written independently; the distance as
  sets of points too (1.7033 / 2.1075 / 3.6942 quadratics per cubic here,
  1.7033 / 2.1074 / 3.6934 in the research run). The kink numbers were redone
  by a second, independent implementation (`kinks.py`); see `kinks.md`.
- **Theorem 3 is wider than the research notes had it.** They had the 4/27
  bound for one kept tangent. The note also proves, the same way, that no
  weights at all get below sqrt(3)/36, so the mid-point spline is the best
  spline with implied on-curve points under this measure.
- **Other converters are described, not measured.** The research run had
  numbers for FontForge, Skia and cubic2quad from reimplementations. They are
  left out: a number about someone else's software should come from that
  software. FontForge and cubic2quad are described from their own
  documentation; kurbo and lyon from their source.
- **No patents**, by Behdad's rule: none searched, quoted or cited.
- **Not adopted, and said so in the note:** the best weights (0.6% fewer
  quadratics at the default tolerance, and every curve of four or more would
  move); a closed form for the single quadratic at equal parameters (same
  decisions); forbidding the single quadratic for masters (the kinks it makes
  are small beside rounding).
- **Left out of the note:** a closed-form bound on the distance for splines
  (research only, unreviewed); splitting at other than equal parameters;
  FontForge, Skia and cubic2quad numbers; the effect of compiling with fused
  multiply-add.

## For Behdad to decide

- The title, the slug (`cu2qu`), and the wording of the byline and of "About
  this note".
