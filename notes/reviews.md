# Adversarial reviews

Behdad asked for two agents to review the note adversarially once it was
written. This is what they found and what was done about each finding. The
note was already published when they reported; the corrections went out as
two later versions of the page.

## First review: the mathematics (2026-11-23)

A referee agent re-derived every lemma, theorem, proof step, table entry and
the Appendix A listing, by hand, symbolically (sympy), in exact rationals and
numerically, without using this repository's own checks. Its verdict: the
lemmas, theorems, proofs and the listing are right; two errors of substance
in what the text claimed around them, one wrong caption, several
overstatements.

### Errors

1. **The 4/27 bound was claimed for "keeping the tangent"; it is proved for
   keeping the derivative.** Theorem 3's second claim is for s_0 = 1 (the
   first off-curve point is the first candidate; the first handle is 3/(2n)
   of the cubic's). A spline that only keeps the tangent's direction may put
   its off-curve point anywhere on the tangent line, and for a particular
   cubic that can be closer: for (0,0) (1,1) (1,1) (2,0), two quadratics
   reach 0.03337 against (4/27)|d3|/8 = 0.03704, and the single quadratic
   0.25 against 0.296. The note contradicted itself: section 5 rightly used
   sqrt(3)/36, not 4/27, for the single quadratic, which keeps both
   tangents. "cu2qu's weights are the best for two and three" needed the
   same qualifier.
   *Done:* the theorem is stated for s_0 = 1 or s_{n-1} = -1, and the text
   around it (abstract, introduction, section 4, Table 1's caption, Figure
   3) says "starts or ends with the cubic's derivative". A paragraph in
   section 4 says what remains true for tangents: when the cubic's handle is
   across its third difference, sliding along the tangent does not change
   s_0, so a rule that keeps the tangent cannot promise less than 4/27; and
   the candidate is the one point of the tangent line that is right for a
   cubic that is a quadratic. `check_math.py` reproduces the counterexample
   and the arch for which the bound stands.
2. **"Best weights" for seven and nine quadratics were not the best.** The
   table and Figure 3 used one block of five and pairs (0.151879) for every
   odd n from five. The best weights are a linear program's: 0.149384 for
   seven, 0.148243 for nine, 0.148180 for eleven, tending to 4/27.
   *Done:* `best_weights.py` solves the linear program (with exact maxima of
   each piece's cubic as cutting planes); its factors up to fifteen are in
   `cu2qu_math.py` and the figure script; Table 1 (2,296.10 and 4,917.60),
   Table 2's second row at a tenth of a unit (4.306 to 4.291) and the text
   (10.3% to 10.6%) follow.
3. **Figure 2's caption:** the implied on-curve point is on the cubic at the
   middle of every even number of quadratics, not only for two. *Done.*

### Overstatements and gaps, all taken

- "L grows with n (Table 1)": a table is not a proof. The note now gives the
  two growth factors of Phi and that both are below (1 + 1/n)^3.
- The second induction of Theorem 3 said "as before"; its thresholds differ.
  It is now written out.
- The recursive test "used to reject" a tie: not always (it accepted ties at
  54). The second review went further; see there.
- "Nothing like this holds for the distance as sets": now with the referee's
  example of two masters whose quadratics lie on their cubics and whose
  instance halfway is 0.07 off.
- "About a quarter fewer for a renderer": that was against cu2qu; restated
  against HarfBuzz, and again, more exactly, after the second review.
- "The test is not needed": only for the spline. Corrected in the abstract.
- "Every formula and number is checked by a script": toned down, and after
  the second review made more nearly true.

### Wording, all taken

"Equal parameters" used in two senses (pieces of equal parameter length; the
same parameter on both curves); the range of e^2 in section 7 (zero when e
changes sign); "smooth exactly when" the cross product is zero (only when);
the instance bound needs no negative coefficient; the introduction's items
for sections 4 and 5 lacked "at equal parameters"; "ordinary input" (finite);
the exact condition behind "add up to 4/3"; Lemma 2 needs pieces of equal
parameter length; distances are Euclidean; ties are "to within rounding";
symbols used twice (a, b; e; "weights"); "a point at an extremum stays an
extremum" (its tangent stays level or upright); "behind a control point"
(behind the end point).

## Second review: claims, scripts and the page (2026-11-23)

A second agent, with four helpers, checked everything but the proofs: every
number against a rerun, every claim about other software against that
software, every cited source against its text, and the page in a browser at
desktop and phone widths, by mouse and by touch, and in print. Its summary:
the numbers from the three fonts hold almost without exception; what was
wrong was elsewhere. The related-work section under-credited its sources,
several sentences about kurbo, HarfBuzz and the masters said more than the
evidence, and the page had layout and touch bugs.

Every finding below was checked again here before anything was changed, with
a script that is now in `scripts/` where one was needed.

### Wrong

1. **"kurbo's port overflowed its stack in the same way" (section 5).**
   Upstream kurbo returns from the two inputs the paragraph listed. It
   overflows on others: a spline at its limit, and single quadratics that
   tie in its own arithmetic. And fontTools fails in the search over n too,
   which the note said only of kurbo. So "the test could reject a tie" was
   one of three outcomes: accept, reject, not return.
   *Done:* `check_endless_kurbo.py` (new) and `check_endless.py` (two more
   inputs) show which input fails where; `ties.py` (new) counts the three
   outcomes on random integer cubics exactly at a limit (54: all 3,000
   accepted; 128: 2,145 of 3,000 rejected; 843.75: 2,290 of 4,200 rejected
   and 61 not returned). Section 5 and the "Ties" remark of section 3 say
   so, with the inputs.
2. **"773,396 groups of compatible masters from 24 variable-font sources,
   of which 15 differ"** was the pull request's number, for a sample that
   also had eight other fonts in it.
   *Done:* `compare_revs_masters.py` (new) runs both revisions on every
   group of the 24: 3,589,652 groups, 142 differ, all in the two CJK fonts,
   all four to three. The note quotes that.
3. **Figure 2's caption** (also the first review's). *Done.*
4. **"cu2qu's weights are the best there are for two and three"** without
   "of those that start and end with the derivative". *Done* with the first
   review.
5. **Table 1's best weights for seven and nine** (also the first review's).
   *Done.*
6. **Table 2, HarfBuzz at half a unit: 2.586 should be 2.585.** The mean is
   2.585454; `check_page.py` rounded a number that `measure.py` had already
   rounded. *Done:* `measure.py` prints six decimals, the page says 2.585.
7. **Reference [8]'s address is dead** (it answers 200 with a page that says
   it cannot be displayed). *Done:* the archived copy is linked, and
   `check_page.py --links` now reads what a link answers with.
8. **"24.4% become one quadratic ... and 32.6% are in fact that close"**
   compared two different quadratics: HarfBuzz's, which borrows a direction
   where a handle has no length, and fontTools', which has none there.
   *Done:* `measure.py` measures HarfBuzz's own: 33.3%. (The reviewer's
   33.6% also counts quadratics for which HarfBuzz's bound gives up, 211 of
   them with handles out of range; 33.3% is of those it considers.)
9. **cubic2quad "with two to eight segments":** its README says so; its code
   tries from one, after cutting at inflection points. The first review's
   one finding not taken; *done* now, from the code.
10. **"24 variable-font sources":** three are fonts as shipped, 21 are design
    spaces of 19 Noto families. *Done;* `kinks_inputs.py` (new) lists them
    with versions in `data/kinks-inputs.txt`.
11. **"That takes a feature smaller than the tolerance"** (a tangent turned
    around by the single quadratic). fontTools returns
    [(0,0), (-2,0), (8,-2)] for (0,0) (2,0) (-2,0) (8,-2) at a tolerance of
    1: a handle twice the tolerance, turned around.
    *Done:* the right limit is 2 sqrt(3) = 3.46 times the tolerance, because
    a quadratic that passes at equal parameters has |e1| and |e2| within
    that; `check_math.py` checks the bound and the example.
12. **"Where the handles of two masters cancel"** (the 69 joints): at 19 of
    them nothing cancels, the handle has no length in any master.
    *Done:* `kinks.py` now says how the 69 split (50 and 19), and the note
    with it.

### Unsupported or overstated

1. **Related work and what is new.** The sources hold more than the note
   credited. Truong, Yuksel and Seiler measure at equal parameters and give
   1/54 of the third difference for the equal pair, which is Phi(2)/8, and
   suggest cutting the cubic by that error. Colomitchi has the measure, the
   two candidates and that they are half the third difference apart, the
   4/27, the mid-point quadratic's error and the cube law, and computes his
   cuts from it. FontForge's documentation already says that with an even
   division the on-curve points fall midway and can be left out. And the
   research notes' own list of what to cite had not reached the section.
   *Done:* section 9 is rewritten from the sources, read again: Colomitchi's
   archived page, Truong et al.'s paper, FontForge's and cubic2quad's text
   and code, and Biły's 2014 paper, which was fetched and read for this
   (mid-point quadratic, error, cube law, and the number of parts as a cube
   root). Watkins and Worsey 1988, Lachance 1988 and Eck 1993 are cited for
   degree reduction with the note saying that they were not read. The last
   paragraph says what is added here, narrowly: the error of the spline for
   every n and every choice of weights, its largest value for the ramp from
   three on, Theorem 3 for every n, and Theorem 4. `check_math.py` checks the
   constants (1/54, 432, 1/32) and that FontForge's observation is exact for
   the halfway points and only near for the meeting points of tangents
   (off by 4 units on an arch 320 wide).
2. **"Every formula and number is checked by a script".** Not checked:
   section 6's numbers, the font versions and "1.1%", the numbers quoted
   from pull requests, the kurbo cubic, and "the closed form decides every
   cubic as the recursive test does".
   *Done:* `check_page.py` now reads `kinks.txt`,
   `kinks-within-one-degree.txt`, `ties.txt`, `endless.txt`,
   `endless-kurbo.txt`, `compare-kurbo.txt`, `perf-kurbo.txt`,
   `compare-fonttools-masters.txt` and `cubics.fonts.txt`; `compare_revs.py`
   compares the closed form of the single quadratic's error with what
   fontTools returns (no cubic of 77,533 differs, at three tolerances). The
   numbers that were quoted from pull requests are replaced by measurements
   made here, except HarfBuzz's "14 to 19%", which now says whose it is.
3. **The population behind 616 / 70 / 26.** "Exactly smooth in every
   master" leaves out nearly every joint off the axes; 4,921 of the 4,939
   are in the two CJK fonts; "the largest change is 4.4 degrees" is at the
   three sampled places.
   *Done:* the paragraph says how a joint was counted and where the joints
   are, gives the counts for joints smooth to within a degree (34,785 off
   the axes beside a single quadratic: 4,951 / 547 / 33 changed by
   conversion, 28,885 / 22,525 / 18,960 by rounding; `kinks.py
   --smooth-within 1`, output now in `data/`), and 4.6 degrees for the worst
   place between that joint's masters. The reviewer's count of joints that
   designers marked smooth (18,095 off the axes, 126 of them exactly smooth)
   was not measured again: the arrays do not carry the marks.
4. **"The conversion changes nothing ... except at 69".** At 29 of the 69
   it changes the angle by more than a tenth of a degree, up to a half turn;
   two more joints change in a master itself. *Done:* said, with why.
5. **kurbo "from 2,426 to 735"** was the pull request's benchmark, which
   alternates two accuracies. *Done:* `perf_kurbo.py` (new) counts each:
   2,228 to 714 at a tolerance of one unit, 3,180 to 756 at half a unit.
   (The reviewer, with the pull request's binaries, got 2,164 and 3,063
   before; the builds differ.)
6. **"Bit for bit"** needed "on these fonts": within rounding of a limit
   the two can differ by a quadratic, since kurbo takes the third
   difference's length as sqrt(x*x + y*y) and Python as hypot. *Done:*
   `compare_kurbo.py` (new) shows both: no bit differs on the 77,533 cubics
   at three tolerances, and one cubic at 128 where kurbo gives 3 and
   fontTools 4.
7. **"About a quarter fewer" for a renderer.** From HarfBuzz as it is,
   accepting its splines by measured distance saves 10.8%, 14.9% and 19.6%.
   *Done:* `measure.py` computes exactly that; the note says 11%, 15%, 20%.
8. **"Most of the difference is in the single quadratic":** 69% at a
   tolerance of one unit, 52% at a half, 22% at a tenth. *Done:* said.
9. **"Add up to about 4/3 ... however that is split":** for very unequal
   handles the sum is not enough. *Done:* the sentence gives the correction
   and says that it grows with the square of the difference.
10. **Abstract: "the test is not needed", "ordinary input".** *Done* with
    the first review.
11. **"So is every instance":** only with no negative coefficient. *Done*
    with the first review.
12. **Table 2's "HarfBuzz" row** is its converter given every cubic;
    HarfBuzz draws a nearly straight cubic as a line first. *Done:* the row
    is named "HarfBuzz's converter", and `measure.py` counts the lines
    (2.1% of the cubics at half a unit; 2.585 becomes 2.526).

### Bugs in scripts and in the page

1. **Theorem 4's formula was cut** (945 pixels in a column of 900), and in
   print four more and the listing. *Done:* the wide formulas are set on
   two or three lines; in print the listings are smaller and nothing is
   clipped. `check_layout.py` (new) fails if a display formula, table or
   figure is wider than its box at print or desktop width.
2. **Touch dragging failed:** the browser scrolled the page and cancelled
   the pointer. *Done:* the control points take `touch-action: none` with a
   `display` that makes Chrome honour it, and a touch that starts on one is
   kept from scrolling. `test_touch.js` (new) drags with real touches in a
   Chrome that emulates a phone. Run on the page as it was, it fails eight
   times; on the page as it is, it passes. Still not tried on a phone.
3. **The page was wider than a phone:** one inline formula, 500 pixels,
   could not wrap. *Done:* it is a display formula now; `check_layout.py`
   checks widths of 320, 390 and 480.
4. **Without scripts the text is raw TeX, with no notice.** *Done:* a
   notice. The formulas are still TeX without scripts.
5. **`check_page.py`:** rounded twice; passed a dead link; nothing compared
   Figure 3's percentile marks with the data. *Done,* all three.
6. **`measure.py`:** the check that the bound is never below the measured
   distance left out the cubics where HarfBuzz borrows a direction. *Done:*
   all 73,924 cubics the bound applies to are in it.
7. **Figure 5 said "the bound does not apply" where it does** (a handle of
   0, and 2.00 by rounding). *Done:* the bound is a function of the handle
   fractions in both `cu2qu_math.py` and the figure script, and the figure
   calls it with the fractions it was given.
8. **Smaller things on the page.** Figure 3 could not be read by touch or
   keyboard, and its readout came and went: it has a slider now, the
   readout stays, and its room is kept. Any mouse button started a drag,
   and a press whose release was lost left the point following the pointer:
   only the main button drags, and a mouse that moves with no button down
   lets go. No figure could be reset: each has a button. Figure 1's slider
   could not come back to a tolerance of 10: its scale now has 10 in the
   middle. In Figure 5 a control point left the drawing past 1.22: the
   drawing has room up to 1.33 and the sliders stop there. The controls
   printed: not any more. The figure script was cached for hours: it is
   linked with its version, as the stylesheet is.
9. **Kinks cannot be rerun from the repository.** The masters are still not
   in it (272 MB); `data/kinks-inputs.txt` now says what they were made
   from, with package versions and the commit of the sources.

### Unclear or minor

- "ended" for the 2015 quotation, which a question follows; and the sketch
  divided by n where cu2qu divides by n - 1. *Done:* "left one thing open",
  and a sentence under the formula.
- "Smooth exactly when": *done* with the first review.
- Reference [6] had a made-up title. *Done:* the pull request's title.
- "The repair in both is": the pull requests are open. *Done:* "proposed".
- "Up to 100": that is fontTools. *Done:* said.
- Figure 1 never shows one quadratic. The figure still draws the spline,
  from two up; its readout now says where one quadratic would do, and the
  caption says that it does. A readout could read "= 128.0: cu2qu 4":
  *done,* it shows as many digits as put the number on the right side.
- The contents list lacked "About this note". *Done.*
- Not changed: pull request 6882's text says 24.1% where the merged code
  gives 24.4% (the text is of an earlier commit; the note uses 24.4%);
  Table 3's numbers differ from the fontTools pull request's own (another
  workload). kurbo issue 607 is not mentioned in the note.

### Found while fixing

- A browser breaks a line between an inline formula and the comma after it.
  The page now keeps a formula together with the punctuation that touches
  it.
