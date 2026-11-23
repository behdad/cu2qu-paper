# How the closed form came about (HarfBuzz)

Notes from the HarfBuzz work, 2026-11-22 and 23, as written then.

PR https://github.com/harfbuzz/harfbuzz/pull/6879, merged 2026-11-23 as 0c72f9596 (branch `gpu-cu2qu-spline`,
opened 2026-11-22). Started as "let's do cu2qu" after the timeout fixes
(`project-gpu-cu2qu-unencodable`, `project-raster-nan-flatten`).

**What the old code did:** `src/hb-gpu-cu2qu.hh` applied fonttools' *single*
quadratic step (control point = intersection of the end tangents)
recursively by halving, depth 10. Three costs:
- the error test (`hb_gpu_cubic_farthest_fit_inside`) is exact but compares
  at equal t; a tangent-intersection quadratic runs at another speed, so it
  over-reports (3x for a quarter circle) -> most cubics one level deeper;
- halving only gives 2/4/8 pieces;
- a control point on its end point (c2 == c3 or c0 == c1) gives no tangent,
  the fit fails down that whole spine to the depth cap (NotoSansSC 3,473
  pieces, NotoSansCJKtc-VF 17,988). The "parallel tangents" cap hits were
  tiny pieces of the same spines (absolute 1e-12 threshold), not
  inflections, as I first guessed.

**Behdad's steer:** he approved "zero-length tangent takes its direction from
the other control point" as "the correct fix", said "you might need to also
fix the error function" and "the error function indeed measures t for t ...
without sampling though". The zero-tangent fix alone removes the cap hits
but saves only 0.05-0.3% of blob bytes, because the equal-t error at that end
is first-order in the piece size.

**Measured options** (curves per cubic on NotoSansSC / SourceCodePro /
STIXTwoMath / NotoSansCJKtc-VF, at upem scale, tolerance 0.5):
- base 4.55 / 4.22 / 5.79 / 5.59
- zero-tangent only 4.53 / 4.22 / 5.68 / 5.51
- fonttools spline (shipped in the PR) 2.62 / 2.59 / 2.99 / 3.04; blob bytes
  -14..-24%, encode-all instructions -18..-31%; tangents exact (0.000 deg at
  3.2M joints and 3.3M cubic ends)
- closed-form kurbo/lyon split (n from third difference, control point
  (3(p1+p2)-p0-p3)/4): 1.82 / 1.92 / 2.07 / 2.07, bytes -20..-32%,
  instructions -41..-55%, C1 inside a cubic but 0.6-1.2 deg off the cubic's
  tangent at its ends on average (>5 deg at 1-2%). Not taken; offered in the
  PR as the alternative.

**Design facts worth keeping:**
- Trial n starts at the provable lower bound: at equal t no quadratic is
  closer to a cubic than |third difference| / 32 (Chebyshev), so n^3 >=
  |d3| / (32 tol). NaN-safe loop; n > 64 -> emit control polygon (that
  replaced #6875's split-point check; its tests still pass).
- Real fonts need n <= 10; random cubics with control points in +-8000 up
  to 34.
- The gain partly comes from using the tolerance: old output sat well
  inside 0.5.
- At 30 units/em (fuzzer scale) 5,133 blobs grow by <= 160 bytes: an
  S-curve the old code did with two tangent-intersection halves now takes 3
  spline pieces.
- No CPU-side decoder of gpu blobs exists in the tests, so fidelity was
  checked at the conversion with temporary counters (sampled equal-t
  deviation max 0.4998), plus byte-identity for all non-CFF fonts.

**Tools** (`memory/gpu-tools/`): `cu2qu_stats.py` (counters for the old
code), `cu2qu_stats2.py` (for the spline code), `c2q_fix_v1.py` /
`c2q_fix_v3.py` / `c2q_fix_v5.py` (the three prototypes as patch scripts),
`c2q_variant.sh`, `c2q_v3.sh`, `c2q_time.sh`, `run_c2q.sh`, `c2q_probe.cc`
(units and blob bytes for sample cubics), `gdiff_c2q.py` + `run_gdiff2.sh`
(gpu_dump marks faces with CFF/CFF2), `gdiff_larger.py`, `fuzz_gpu.sh`.

**Paper idea (Behdad, 2026-11-22):** "we should write a paper on it some time" (he wrote the cu2qu spline algorithm). Material from this session that would feed it: the three reasons recursive tangent-intersection wastes curves; the |d3|/32 Chebyshev lower bound on n; the measured comparison with the closed-form split (curve counts, tangent error at cubic ends); my unverified arithmetic that a one-sided end piece has equal-t error (4/27)|d3| against sqrt(3)/36 |d3| for the midpoint control point. No commitment or timeline.

**Direct n shipped (2026-11-23):** HarfBuzz PR https://github.com/harfbuzz/harfbuzz/pull/6880 merged as e6ce4cc7b (n from the closed form, no trial, no per-piece test; 3 exact-tie glyphs of 137,440 CFF pairs got one quadratic fewer; -13.6..-18.8% encode instructions). HB's product form `d3 * Phi(n) <= tol * n^3` accepts exact ties at 54 and 128, same as the fonttools table form. fonttools side: see `project-fonttools-cu2qu-compute-n`.

**Cleaner closed form (found while writing the fonttools docstring):** every interior piece k (0 < k < n-1) has the same shape, `phi_k(v) = (n-1-2k)/(4(n-1)) + c v - v^3` with `v = u - 1/2`, `c = (n-3)/(4(n-1))`: a fixed odd cubic plus a constant ramping linearly from +c to -c. Worst interior piece is k = 1: `c + 2 (c/3)^1.5` at `v = sqrt(c/3)`. First piece `u^2 (1+g) - u^3`, max `(4/27)(1+g)^3` at `u = 2(1+g)/3`. The first piece is the worst up to n = 54, the second from n = 55. Why it all works: the two points `cubic_approx_control` goes between are d/2 apart, and the halfway points of consecutive pieces have the cubic's split point as their midpoint.

**Exact n for the spline (worked out 2026-11-23 at Behdad's "can you work out a lower bound on n?"), verified by `memory/gpu-tools/spline_bound.py`:** for n >= 2 the difference between quadratic k and its piece at equal u is `(d3 / n^3) * phi_k(u)`, a multiple of the one vector d3, with `phi_k(u) = g_k - u (2 g_k + k/(n-1)) + u^2 (1 + g_k + g_{k+1} + k/(n-1)) - u^3`, `g_k = (n-2k)/(4(n-1))` for 0<k<n and g_0 = g_n = 0 (g_k d3/n^3 is the offset of implied on-curve point k from the cubic). So n fits iff `|d3| * Phi(n) <= tol * n^3`, Phi(n) = max |phi_k|: Phi(2) = 4/27, Phi(3) = 27/128, increasing to 1/4 + sqrt(3)/36 = 0.29811; closed form = max of the first piece `(4/27)(1+g_1)^3` and the second `(n-3)/(4(n-1)) + 2((1/4 - 1/(2(n-1)))/3)^1.5` (or g_1). Thresholds n^3/Phi(n): 54, 128, 272.05, 503.86, 843.75, ... Lower bound n >= cbrt(4|d3|/(27 tol)); n >= cbrt(0.29811 |d3|/tol) always fits. Checked: identity to 8e-15 on 2,000 random cubics; fontTools trial n equals the formula on all 196,545 random cubics that needed more than one quadratic. Not yet used in HB (PR 6879 starts trials at the looser |d3|/32 Chebyshev bound); using it would remove the trial loop and the per-piece test with identical output. Venue for a write-up: he said "it can totally just be on my website".
