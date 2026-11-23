# The fontTools pull request

Notes as written on 2026-11-23.

PR https://github.com/fonttools/fonttools/pull/4452, branch `cu2qu-compute-n`
pushed to fonttools/fonttools (Behdad: "PR both fonttools & HB", "~/fonttools
is yours"). Two commits: 41306dff0 "[cu2qu] Compute the number of quadratics
instead of trying", 6f811387b "[cu2qu] Do not recurse without end on a NaN
tolerance". HarfBuzz half: `project-gpu-cu2qu-spline` (PR 6880, merged).

**Design:** `spline_error_factor(n)` is a plain untyped `def` (no pow, no typed
C arithmetic) so the table `SPLINE_MAX_D3[n] = n^3 / error(n)` is the same bits
compiled or pure and on every platform (54.0 and 128.0 come out exact);
`cubic_spline_n` (cfunc) returns the least n with `abs(d3) <= tol *
SPLINE_MAX_D3[n]`, 0 if none up to MAX_N; `cubic_spline(cubic, n)` builds with
the same point arithmetic as before (do not "optimize" it: low-bit changes
flip otRound at .5, they have tests for that). The single quadratic is still
tried first with the recursive test. `curves_to_quadratic` takes the max n
(limits grow with n).

**Output differences (I first built a byte-identical conservative variant;
he said output may change for the better, see
`feedback-output-may-change-for-better`):** only cubics within rounding of a
limit. Real fonts: 10 of 155,066 conversions, all exact ties (e.g. |d3| = 128
tol took 4 quadratics, now 3). Synthetic at-the-limit class: 7,309 of 40,000;
checked in 60-digit arithmetic, new is right except within 1.1e-15 of the
limit, old recursive test was off by up to 1e-9 (large coordinates).
Compiled == pure on all 515,066 inputs, before and after.

**NaN:** old compiled module segfaults (pure: RecursionError) on NaN
coordinates and NaN tolerance: 64 of 280 odd inputs. Commit 1 fixes
coordinates (d3 NaN -> ApproxNotFoundError), commit 2 the tolerance
(`if not abs(mid) <= tolerance` in cubic_farthest_fit_inside).

**Speed (us per call, compiled / pure):** font cubics tol 1: 1.53 -> 1.23 /
5.7 -> 3.5; tol 0.5: 2.09 -> 1.48 / 8.1 -> 4.3; cu2qu.benchmark
curve_to_quadratic 45.5 -> 14.4, curves_to_quadratic 98.1 -> 48.9.

**Third commit (fe40caaf1), after he asked "for multi-curve API doesn't use
the new skipping bound eh?":** `cubic_approx_quadratic` gives up before
intersecting tangents when `d3 > 21 * tolerance` (`QUADRATIC_MAX_D3`): any
quadratic with the cubic's end points is at least sqrt(3)/36 |d3| away at
equal parameters (36/sqrt(3) = 20.78), so results are unchanged. A first
version that hoisted the check into both callers with a `d3s` list made the
compiled multi-curve path slightly slower; the check inside the shared
function gives -6% (single) / -2.4% (multi) compiled, -9..-11% / -2.6% pure,
in instructions. Wall-clock benchmarks were useless that day (machine loaded
by subagents, 2-4x swings): use `perf stat -e instructions:u`
(`ft_perf.py` / `ft_perf.sh`).

**kurbo follow-up:** done as a branch, see `project-kurbo-cu2qu-port`.

**Verification setup:** fonttools CI = tox lint (black==26.3.1 pinned, mypy
only on plistlib), py-cov, py-cov-noextra, py-cy (Cython), pypy3. Locally: venv
with Cython 3.3.1 + requirements.txt for python3.14 at
`~/.cache/cu2qu-work/venv` (deleted 2026-11-23 with the build trees; recreate with `ft_venv.sh`, pointing its W at a directory under ~/.cache), exports of old/new built in place with
`FONTTOOLS_WITH_CYTHON=1 python setup.py build_ext --inplace` (never build in
place in ~/fonttools: a stale .so shadows the .py). pytest there runs
doctests too (`--doctest-modules --pyargs`, testpaths Tests + fontTools).
Scripts in `memory/gpu-tools/ft-cu2qu/`: `ft_equiv.py` (fixed-seed dump of
results, SKIP_ODD / TIE_INPUTS env), `ft_tie_exact.py` (60-digit check of
each difference), `ft_odd.py` (one process per odd input, crash = outcome),
`ft_bench.py`, `ft_build_cy.sh`, `ft_dumps.sh`, `ft_fulltests.sh`,
`ft_verify_all.sh`, `ft_old_tests.sh` (new tests against old code).

**Two follow-up agents (2026-11-23, his "do those in agents and update PRs
with findings"):**

*Kinks from the single quadratic in masters: leave it.* His requirement:
"condition for kink-free interpolation is either equal angle or equal ratio of
handles", and exact end tangents are "nonnegotiable for fonts". n >= 2 is
linear in the control points, new on-curve points have ratio exactly 1, end
handles are 3/(2n) of the cubic's: proved and checked on real masters. n = 1
(tangent intersection) gives handles 1/l times the cubic's with l per master,
so equal ratio is not kept: tan(kink) = (R - 1) sin(theta) / (2 sqrt(R) +
(1 + R) cos(theta)). Measured on 24 real sources (Cantarell, Noto Sans/Serif
CJK CFF2 masters, 21 noto-source design spaces; 3.59M cubics, 448k single
quadratics, 194,218 smooth joints next to one): angle changed by more than
0.1 / 0.5 / 1 degree at 645 / 77 / 27 joints, against 3,982 / 3,387 / 3,062
from integer rounding; worst real case 4.6 degrees (4.5-8 unit handles, 34
degree turn between masters). Forbidding single quadratics for masters costs
+2.4% points (8% Noto Naskh, ~25% for TrueType-derived masters). No code
changed. By-product: kurbo branch == fonttools branch bit-identical on
773,396 real master groups; fonttools main vs branch differ on 15 (4 -> 3).
Scripts: `memory/gpu-tools/cu2qu-research/kinks/`.

*The fit test recurses without end on finite input, in released fontTools.*
`curve_to_quadratic([(0,0),(0,0),(0,0),(25,0)], 100/27)` and
`curve_to_quadratic([(-45,-3),(-30,-2),(0,5),(0,9)], 20/3)`: RecursionError
pure, segfault compiled (I reproduced on 4.60.1 and the branch). Mechanism
(agent's proof): piece (A, A+ulp, A, A) with tolerance A, A's significand
3 mod 8: mid rounds to A, deriv3 = -ulp/8, first half is the same piece.
Needs the error curve to stop at its farthest point (three equal points, or
farthest point on an axis), which my earlier generators skipped. Fix: depth
limit 64 -> True, as kurbo (commit 4 on PR 4452; kurbo got the inputs as
tests, commit 5f6ba66). Tests that end on their own stay under 40 halvings
(70,000 at-the-limit inputs, mine). Closed form via the cubic equation was
evaluated and rejected: same decisions on 77,533 font cubics, same work,
more code, overflow guards needed.

*FMA side finding, not acted on:* built with `gcc -mfma` the compiled cu2qu
has 134 fused multiply-adds and differs from pure Python on non-integer
input (2,305 of 120,000 random cubics in coordinates, 11,720 of 210,000 tie
decisions); integer input identical; default x86-64 identical. Whether arm64
wheels are built with contraction is unchecked; `-ffp-contract=off` would be
the cure. Told Behdad; no issue filed.
