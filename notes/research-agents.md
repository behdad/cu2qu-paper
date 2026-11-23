# Research behind the note

Notes taken on 2026-11-23 from three research runs (prior art and other converters; a better
error measure; better weights), as they were written then. Where the note and these differ, the
note has what was checked again by the scripts in `scripts/`; the rest here was not re-verified.
Patents are left out on purpose: none is searched for, quoted or cited, here or in the note.

Behdad wants a paper / technical note on cu2qu ("we should write a paper on it
some time"; venue "JCGT sounds good" or "it can totally just be on my
website"; motivations "Slug algorithm for one", "TrueType fonts for other",
"with truetype the algorithm has the nice property of resulting in the implied
mid oncurve point"). On 2026-11-23 he had me "spawn three research agents, one
for optimal / great weights, two for better error measure (analytical, no
sampling), and three to study kurbo and other legit approaches and their
results". Math background: `project-gpu-cu2qu-spline`.

## Survey agent (finished; its tags: read in source / measured / recalled)

**Prior art for "n in closed form":** only for the mid-point scheme (s = 0):
kurbo `to_quads` (uniform closed form since 2020-01, Raph Levien), lyon >= 1.0
(`num_quadratics`), Vello classic, Impeller, Skia Graphite sparse strips; in
print Biły 2014 (Theorems 3, 6, sec 5.4.1); Colomitchi 2006 has the constant
sqrt(3)/36 but iterates. **Nothing found for cu2qu's ramp at general n**:
fontTools, kurbo `approx_spline`, fontc, FontForge, cubic2quad all search n
by trial. Only n = 2 is known: Truong-Yuksel-Seiler, HPG 2020, build exactly
cu2qu's n = 2 spline (gamma = 1/2) with max error |a3|/54 (eq. 29); the 4/27
is Colomitchi's eq. (5). Caveat: Eck 1993, Watkins-Worsey 1988, Lachance 1988,
Bogacki et al. 1995, Sederberg 10.6 were not read (paywalled).

**My "alternating" LP optimum is prior art:** Truong pairs on n/2 equal
pieces: implied midpoints, tangents kept, (4/27)|d3|/n^3 for every even n
(limits 54, 432, 1458, 3456 vs cu2qu 54, 272.05, 843.75, 1930.24). n = 4 is
Helen Triolo's "fixed MidPoint" routine published by Groleau 2002.

**Chebyshev:** |d3|/32 is standard economization (Watkins & Worsey 1988,
Lachance 1988, Eck 1993; Ahn 2003 "well-known"). With end points interpolated
the best is sqrt(3)/36 |d3| = |d3|/(32 cos^3(pi/6)), constrained Chebyshev
t(t-1/2)(t-1): the mid-point quadratic is optimal there.

**Measured on the 77,533 font cubics (quads per cubic at tol 1.0 / 0.5 / 0.1;
two-sided Hausdorff, checked against brute force to 2.5e-5):**
- cu2qu = kurbo `approx_spline` = HarfBuzz: 2.336 / 2.838 / 4.802; true
  Hausdorff averages 0.27 / 0.34 / 0.45 of the tolerance.
- kurbo `to_quads` = lyon 1.0: 1.688 / 1.984 / 3.037, but end tangent off by
  more than 1 degree on 57.4% of cubics; bounding box changed on 31.8% of
  monotone cubics with an axis-aligned end tangent (up to 0.77 units).
- cubic2quad (svg2ttf; tangent intersection per piece, sampled test):
  1.640 / 1.860 / 2.544 quads, 2.253 / 2.695 / 4.040 TrueType points
  (on-curves not implied), tangents kept at every joint.
- FontForge (own tolerance 1.5 with rounding): 1.607 quads, 1.593 points.
- Skia convertCubicToQuads at an equivalent bound: 2.488 / 3.149 / 5.314.
- Oracles (least n by true Hausdorff): cu2qu's construction 1.703 / 2.107 /
  3.694; mid-point 1.417 / 1.698 / 2.611; tangent-intersection off-curves with
  implied on-curves 1.633 / 1.907 / 2.943.
- **The equal-parameter test is the main cost**, worst for the single
  quadratic: cu2qu accepts one quadratic for 4,893 cubics at tol 1.0, while
  the same tangent-intersection quadratic is within 1.0 geometrically for
  35,496.
- Tangent-intersection quadratics are 4th-order accurate geometrically
  (measured 4.00; mid-point 3.00; ramp tends to 3). Theory: Schaback 1989,
  Floater 1995.
- Distribution at tol 1.0 (n = 1/2/3/4/more): cu2qu 6.3 / 65.9 / 17.6 / 8.7 /
  1.6 %; kurbo 38.6 / 54.1 / 7.0 / 0.2 / 0 %.

**Holes in cu2qu's "tangents kept":** a zero-length handle (862 cubics) puts
the first off-curve point on the end point: tangent off > 1 degree on 641,
> 10 on 331. After rounding off-curve points to integers (as glyf stores)
9,790 cubics (12.6%) are off by > 1 degree and 848 exceed tol 1.0 (max 1.45).

**HarfBuzz vs fontTools on the corpus:** counts differ on 3 / 13 / 0 cubics:
the |d3|/tol = 128 ties (closed form 3, trial 4) and 6 single-quadratic ties
at tol 0.5 where cu2qu returns 1 and HarfBuzz 2 (agent's guess: HB's depth cap
of 8 in the recursive test; not verified by me).

**A write-up must cite:** fontTools cu2qu (2015; ramp in commit e686d61,
ancestor Roboto's 2014 convertCurves.py); Truong-Yuksel-Seiler 2020;
Colomitchi 2006; Biły 2014; kurbo and lyon; Watkins & Worsey 1988, Lachance
1988, Eck 1993, Bogacki-Weinstein-Xu 1995; FontForge (PrettyApprox blends two
tangent chains linearly in the piece index); Triolo / Groleau 2002; Alfke
1994; Apple TrueType Reference Manual Appendix E (warns that comparing at the
same t "may not work ... different velocities", prefers perpendicular
distance); Schaback 1989 or Floater 1995; cubic2quad; Lengyel 2017.
(No patents in this list or anywhere in the write-up: see
`feedback-no-patents`.)

Files: `memory/gpu-tools/cu2qu-research/survey/` (Rust harness
`bench/src/main.rs` + `reimpl.rs`, `run_all.sh`, `report_tables.py`,
`analyze_percubic.py`, `out/*.txt`), corpus `cu2qu-research/ft_cubics.json`.
Large data (1.4 GB, regenerable) was in `~/.cache/cu2qu-work/research-survey/`.

## Error-measure agent (finished; proofs are its own, numerically checked, not reviewed)

**Headline (corpus, cu2qu ramp, quads per cubic at tol 1.0 / 0.5 / 0.1):**
equal-parameter today 2.3364 / 2.8380 / 4.8016; closed-form geometric bound
with equal-parameter as fallback 1.7336 / 2.1370 / 3.7178 (-25.8 / -24.7 /
-22.6 %); true least n by measured Hausdorff 1.7033 / 2.1074 / 3.6934. The
bound was never below the measured distance in 2,237,419 (cubic, n) pairs
(corpus + 20 random/adversarial sets) and never picked fewer quadratics than
the truth. Reference implementation ~330 lines of plain Python
(`cu2qu_geo.py`), 5x slower than pure-Python fontTools as written; no C port.

**n = 1 identity (I re-checked it: 6e-14 on 2,000 random cubics).** With q the
tangent intersection, a = q - c0, z = c3 - q, c1 = c0 + l1 a, c2 = c3 - l2 z,
p = l1/2 - 1/3, r = 1/3 - l2/2, e(t) = (1-t) p + t r:
`B(t) = Q(tau(t)) + rho(t) (a - z)`, `tau = t + 3 t (1-t) e(t)`,
`rho = 3 t^2 (1-t)^2 (l1 + l2 - 4/3 + 3 e(t)^2)`. The cubic leaves the parabola
only along its axis direction a - z. So a single quadratic fits when the two
handle ratios ADD UP to 4/3, however unevenly split; the equal-parameter test
wants each to be 2/3. Rule of thumb: distance ~ (3/8) |l1 + l2 - 4/3| h, h the
distance of q from the chord (median ratio 1.000). At tol 1.0, 45.8% of corpus
cubics are within tolerance of their tangent-intersection quadratic; today's
test accepts 6.3%; the closed-form bound accepts 43.7%.

**n >= 2.** Error across the spline's tangent is
(|d3|/n^3) |phi_k(u)| |sin(angle of d3 to the tangent)|, and d3 x B' is linear
in t. First-order term alone is NOT a bound (0.911 of the truth for n = 2;
would pick too few for 833 corpus cubics). Lemma N gives the exact distance
along the normal (reproduces Hausdorff to 1e-5); Theorem 2 bounds it by
T1 + T2 + T3 per piece with validity conditions; fallback to equal-parameter
where they fail (0-1.5% of corpus cubics; most cusps, loops, coincident
points). Formulas in its report and in `closed.py` / `cu2qu_geo.py`.

**Equal-parameter test for n = 1 in closed form:** max of
|3u(1-u)^2 e1 + 3u^2(1-u) e2| is at a root in (0,1) of
`3c u^3 + (5b - 2c) u^2 + (2a - 3b) u - a = 0` with a = |e1|^2,
b = e1.(e2 - e1), c = |e2 - e1|^2 (I re-derived this). So cu2qu could be
entirely closed form: Cardano for n = 1, Phi(n) for n >= 2, no recursion.

**Interpolatable masters: the catch.** Equal-parameter error is convex in the
control points (masters fit => every instance fits); geometric distance is
not. With the n = 1 geometric test per master, the half-way instance was
outside tolerance in 108 of 1,986 pairs at tol 1.0 (worst 3.05x); for n >= 2
only 2 of 20,000 (worst 1.04x). So `curves_to_quadratic` cannot use the n = 1
geometric test per master as is; fine for single-instance rendering (hb-gpu).

**Other weightings with the geometric bound (corpus):** s = 0: 1.4728 /
1.7492 / 2.6387 (GPU case at 0.5: 1.749 vs 2.838 today, -38%); alternating:
1.8635 / 2.3082 / 3.7216.

Files: `memory/gpu-tools/cu2qu-research/error/` (`cu2qu_geo.py`, `closed.py`,
instrument `gt.c`, `final_tables.py`, `n1_identity.py`, `check_*.py`,
`interp_*.py`, tables).

## Weights agent (finished; its tags: proved / verified / derived / measured / conjectured)

**Bottom line:** under cu2qu's own (equal-parameter) measure the weights
problem is solved and better weights buy little (-0.6 / -2.6 / -10.6 %
quadratics at tol 1.0 / 0.5 / 0.1), because n = 2 and n = 3 leave no choice
of weights and only cubics with |d3| > 128 tol (10% at tol 1) can benefit.
The large gains come from the error measure, from the tangent constraint,
and from leaving the d3 family for tangent intersections.

**Equal-parameter, uniform splits (my LP table confirmed):**
- Theorem [proved]: with V(x) = max_u u(1-u)(u-h), h = (1-x)/2 (V(1) = 4/27,
  V(0) = sqrt(3)/36), every weighting has Phi >= V(|s_0|). One kept end
  tangent already costs 4/27 for every n; the second costs nothing more for
  even n.
- Even n [proved, also by dual certificate]: +1, -1, +1, ... is optimal and
  unique; it is n/2 sub-cubics with cu2qu's two-quadratic spline.
- Odd n: n = 3 is stuck at 27/128 (optimum (1, 0, -1)); n = 5 optimum
  (1, a, 0, -a, -1), a = -0.933394180511709 (root of a sextic), Phi =
  0.151879364835661; Phi_odd - 4/27 shrinks ~38x per four pieces added.
  Practical odd n >= 5: pairs (+1, -1) plus one block of five, Phi =
  0.151879. Thresholds of |d3|/tol for "pairs": 54, 128, 432, 823.0, 1458,
  2258.4, 3456.
- Relaxing the tangent (s_0 = x): optimum V(x) for even n; near x = 1,
  V ~ 4/27 - (1-x)/9. Allowing the end tangent to turn by a budget, quads per
  cubic at tol 1.0 / 0.5 / 0.1: exact 2.32 / 2.76 / 4.29; 0.5 deg 2.09 / 2.37
  / 3.26; 1 deg 1.99 / 2.18 / 3.10; unbounded 1.70 / 1.99 / 3.05.

**Tangent-intersection spline "T" (off-curve point k = intersection of the
cubic's tangents at knots k and k+1; the n >= 2 form of cu2qu's n = 1 rule):**
keeps both end tangents and the tangent at every knot. Geometric third-order
constant 2/27 at the ends and sqrt(3)/36 inside, against 4/27 for alternating
and Phi(n) for the ramp. Useless under the equal-parameter measure (its
equal-parameter error is 7.6x to 17x its Hausdorff distance). Breaks down near
cusps and across inflections.

**Non-uniform splits:** the right knots are at equal *affine arc length*
(sigma = integral of |B' x B''|^(1/3) dt): an implied-on-curve spline is a
uniform quadratic B-spline, and its third-order error vanishes exactly then.
Verified: the numerically best spline's on-curve points sit at those
parameters to 0.002. q_M (s = 0) on unequal pieces does NOT work even
geometrically; the tangent intersection does. Fourth-order model with
universal constants c(n); closed-form density E0(t) =
h^4 |12 p p'' - 5 p'^2| / (72 |p|^(7/3) |B'|), p = B' x B''.

**Yardstick (numerical optimum of all off-curve points, by Hausdorff), corpus
quads per cubic at tol 1.0 / 0.5 / 0.1:** cu2qu as implemented 2.336 / 2.838 /
4.802; ramp with a geometric test 1.672 / 2.079 / 3.676; alternating 1.667 /
2.059 / 3.367; T at equal-parameter knots 1.603 / 1.879 / 2.912; T at
closed-form knots 1.592 / 1.839 / 2.716; best found with tangents kept (TAN*)
1.543 / 1.753 / 2.299; s = 0 1.417 / 1.698 / 2.611; best with tangents free
(FREE*) 1.144 / 1.253 / 1.782 (85.6% of cubics fit one quadratic at tol 1).

**Its recommendations:** (A) tangents kept, rigorous equal-parameter bound:
"pairs" (-0.6 / -2.6 / -10.3 %). (B) tangents free (GPU), rigorous: s = 0 with
n = ceil(cbrt(sqrt(3)|d3|/(36 tol))) (-27.7 / -30.1 / -36.8 %). (C) geometric
schemes driven by an estimate at 17 fixed nodes per piece, 1.1x margin, 0
cubics over tolerance on the corpus: balanced T with tangents kept 1.653 /
1.899 / 2.719; balanced T with tangents free 1.229 / 1.385 / 2.212 (-47 / -51
/ -54 %). (C) is an estimate, not a bound, and still loops over n.

Not done by it: rigorous geometric bound for T; inflection cubics in the T
family; integer rounding; interpolation compatibility of T (non-linear in the
control points, like the n = 1 rule).

Files: `memory/gpu-tools/cu2qu-research/weights/` (`q1_*.py` LP and proofs,
`q2_*.py` non-uniform and fourth-order, `q3_*.py` sliding, `schemes.py`,
`tquad.py`, `balanced.py`, `gest.py`, `geo2.py`, `opt.py`, `canon.py`, `chd.c`
instrument, `*.out` tables).

## What the three reports say together (my synthesis, 2026-11-23)

- The closed form for the ramp is new as far as the survey could tell, and
  the ramp is within ~10% of the best weights under the same measure.
- The equal-parameter measure is what costs quadratics (~25%), mostly in the
  single-quadratic test; but it is the measure that makes interpolatable
  masters safe (convex in the control points), so fontTools' multi-master
  path should keep it. Single-instance users (hb-gpu) need not.
- For hb-gpu, where end tangents and implied points do not matter: s = 0
  (kurbo's scheme) with the rigorous bound is -30% at tol 0.5; with the
  geometric bound -38%; balanced tangent-intersection schemes about -50% but
  estimate-driven.
