# The kurbo pull request

Notes as written on 2026-11-23.

Branches in `~/kurbo` (on upstream main a5f97ff), pushed to
`git@github.com:behdad/kurbo.git`, only `kurbo/src/cubicbez.rs`:
`cu2qu-halve` = a488b0f "Halve cubics in the cu2qu port as cu2qu does, and
bound the fit test"; `cu2qu-compute-n` = that plus 0ec0ebe "Compute the number
of quadratics in approx_spline instead of trying". He opened it himself as
https://github.com/linebender/kurbo/pull/625 on 2026-11-23 from a draft
description I gave him (see `feedback-third-party-prs`); all 15 CI checks
passed, including wasm32, macOS, Windows and typos that I could not run. It
ports fonttools PR 4452 (`project-fonttools-cu2qu-compute-n`).

Final shape (after his "a bit of diff with fonttools on chronic cases is
fine", 2026-11-23): kurbo's own `Vec2::hypot` everywhere (libm hypot in the
trial loop cost +33% instructions), cu2qu's halving (`subdivide_2`) for
coordinate parity, and a depth limit of 64 in `fit_inside`. The depth limit is
needed even with computed n: the single-quadratic test overflowed at its own
limit, and cu2qu's halving alone did NOT stop it with `Vec2::hypot` (it did
with libm hypot on 835k inputs, but that is luck, not a guarantee: the fixed
point (A, A+ulp, A, A) can occur with any hypot; fontTools itself has no
guard). Versus the fontTools branch: all coordinates identical on 835,066
inputs; counts differ only on engineered ulp-level ties (60 of 40,000 and
4,794 of 280,000). Instructions per call on font cubics: 2,426 upstream,
2,457 after the first commit, 735 after the second.

**What kurbo has:** `CubicBez::to_quads` (its own: equal split, midpoint
control, n = ceil((|d3|^2 / (432 acc^2))^(1/6)), i.e. sqrt(3)/36 |d3| / n^3,
after Colomitchi) and `approx_spline` / `cubics_to_quadratic_splines`, a port
of cu2qu that tried n = 1..=100 (the multi one up to 101) with `fit_inside`,
kept "matching Python" by hand (div_exact, 1/27 reciprocal).

**Findings from a bit-for-bit differential against fontTools (every float
compared by its bits; inputs with integer and non-integer coordinates, at
limits, groups, 77,533 font cubics):**
- Upstream kurbo overflows the stack (process abort) on finite cubics:
  `(386,-158) (527,-28) (480,-99) (2270,2329)` at accuracy 4 (third
  difference (2025, 2700) = 843.75 x 4, the limit for six quadratics; 61 of
  200 random integer cubics with that third difference crash; also with
  quarters at accuracy 1 and eighths at 0.5), and at the single quadratic's
  own limit, e.g. `(-31,97) (-51,93) (-91,76) (-139,52)` at
  1.987615979999813. fontTools does not (7 quadratics, resp. 1).
  Cause: `fit_inside` tests `mid = (p0 + 3(p1+p2) + p3)/8` but splits with
  `ParamCurve::subdivide`, whose middle is `eval(0.5)` and can round to the
  other side of the distance; the halves then never settle. cu2qu builds the
  halves from the tested mid (`mid -/+ deriv3`), which guarantees they end
  on points known to be inside.
- At exact limits upstream is erratic: |d3| = 128 x accuracy gives 3 or 4
  quadratics depending on the rest of the cubic (26 vs 174 of 200).
- Non-integer coordinates: 15,662 of 120,000 random cubics differ from
  fontTools in the last bit of a coordinate, all with n = 2, 4 or 6 (the
  `subdivide` halving in `split_into_n`): 2% / 96% / 86% of those. Integer
  inputs are bit-identical. Simple example: `(5.9,3.6) (0,9.2) (7,6.6)
  (-16.2,-7.7)` at 1.0, third point x 1.199999999999999 vs
  1.1999999999999993.
- `Vec2::hypot` is `sqrt(x^2 + y^2)`, not libm hypot ("slow library
  function"); Python's `abs(complex)` is libm hypot. For decisions to match
  at ulp-level limits both the third difference and the fit test use
  `f64::hypot` on the branch (39 instructions per call).
- After the branch: kurbo == fontTools branch bit for bit on 555,066 +
  280,000 inputs (the latter: single quadratic within 3 ulps of the
  accuracy); no crash. Instructions per `approx_spline` call on font cubics:
  2,426 upstream -> 887.

**kurbo conventions learned:** CI = fmt, taplo, copyright header check,
clippy via cargo-hack (std, no_std+libm on x86_64-unknown-none, wasm32),
nextest, doc tests, MSRV 1.85, docs on nightly with private items, typos.
Lints deny `!(a <= b)` on floats (`neg_cmp_op_on_partial_ord`): write `let
inside = a <= b; if !inside`. No attributes on `if`. no_std needs
`FloatFuncs` for sqrt/hypot. Stable toolchains here: 1.99.0 (their CI's),
1.85. Build with `CARGO_TARGET_DIR` on the home disk.

**Tools** (`memory/gpu-tools/kurbo-cu2qu/`): `kdiff_py.py` (inputs + Python
side), `main.rs` + `Cargo.{old,new}.toml` (Rust side, old = pristine clone,
new = ~/kurbo), `run.sh`, `analyze.py`, `show.py`, `find_crash.py`
(bisects stack overflows), `tie_probe.py`, `gen_single_ties.py`, `one.py`,
`find_simple.py`, `bench.rs` + `bench.sh` (perf stat instructions),
`kurbo_checks.sh` (local CI), `final_check.sh`.
