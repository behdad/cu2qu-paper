# One quadratic by distance (HarfBuzz)

Notes as written on 2026-11-23.

**PR https://github.com/harfbuzz/harfbuzz/pull/6882** (branch
`gpu-cu2qu-geometric-single`, commit 6c0f36fae, on main 77814c49e), opened at
Behdad's "4. do it" after the research agents' reports
(`project-cu2qu-paper-research`). Earlier the same day: PR 6881 (skip the
single quadratic when d3 > 21 tol, an equal-parameter lower bound; -2% encode
instructions, output-identical) was merged as 77814c49e; 6882 removes that
check again because a single quadratic can now fit past it.

**The bound (I proved it myself; the research agent's tighter "bound 2" is not
used).** q = intersection of the end tangents (with HB's zero-handle
fallback), a = q - c0, z = c3 - q, c1 = c0 + l1 a, c2 = c3 - l2 z:
`B(t) = Q(tau(t)) + rho(t)(a - z)`, tau = Bezier function with control values
0, l1/2, 1 - l2/2, 1; `rho = 3 t^2 (1-t)^2 (l1 + l2 - 4/3 + 3 e(t)^2)`,
e(t) = (1-t)(l1/2 - 1/3) + t(1/3 - l2/2). Identity verified exactly in
rationals on a grid covering its degree. For 0 <= l1, l2 <= 2 tau stays in
[0,1] (convex hull) and is onto, so Hausdorff <= (3/16) max|g2| |a - z| with
g2 = l1 + l2 - 4/3 + 3e^2, e^2 between min (0 if p r <= 0) and max of p^2,
r^2. Exact for a symmetric arc (quarter circle: 0.0607 R). Zero handle:
l1 = 0, l2 = 1 gives (3/64)|2c2 - c0 - c3|.

**Numbers (77,533 font cubics):** accepted as one quadratic at tol 1.0 / 0.5 /
0.1: equal-parameter test 6.3 / 2.4 / 0.3 %, bound 37.2 / 24.1 / 8.6 %, truth
about 48 / 33 %. Quadratics per cubic at 0.5: 2.838 -> 2.586 (-8.9%); at 1.0
2.336 -> 1.998. Only 15 cubics at 0.5 pass the old test and not the bound, so
the recursive test (and its depth cap 8) was dropped outright; a hull-bound
fallback recovered nothing. hb-gpu-encode-all instructions vs main with the
skip: -5.6% NotoSansSC, -1.6% SourceCodePro, -1.3% STIX, -2.3% CJK-VF; blob
bytes -1.6% over 445 CFF faces (87,461 of 235,312 glyphs change). Bound vs
measured (point-to-polyline both ways, 18,000 cubics incl. edge ratios, far
intersections, scales 1e-3..1e6, zero handles): never below; median 1.10 on
font cubics, 90% 3.3.

**Why hb-gpu only:** the geometric distance is not convex in the control
points, so masters that fit do not guarantee instances fit (agent: 108 of
1,986 master pairs out of tolerance half-way). Behdad: "geometric [measure]
is a dead end" for fonts, yes; exact end tangents "nonnegotiable for fonts".
fontTools/kurbo keep the equal-parameter test.

**Zero-length handles: measured, not fixed (my recommendation 2 was wrong).**
A cubic with c1 == c0 starts at zero speed with infinite curvature; under the
equal-parameter test a quadratic can follow only with its off-curve point on
the end point (cu2qu's output: a short straight start) or a leg shorter than
about the tolerance, which integer rounding destroys. In the three fonts, of
860 such ends only 10 sit at a smooth join (neighbour within 1 degree; spline
off by median 4.1, max 5.9 degrees), 23 at 1-5 degrees, 827 at corners where
the direction does not matter (`zero_handles.py`). The survey agent's ">10
degrees on 331 cubics" were all corners. In hb-gpu the new single-quadratic
test gives 281 of the 862 an exact tangent anyway.

**A/B pitfall:** `ninja -C build src/libharfbuzz-gpu.so` is a no-op (not a
real target), so my first A/B compared the new code with itself (0.000%).
Build the base with plain `ninja -C build` and `cmp` the two libraries before
trusting an A/B (`hb_skip_ab.sh` is fixed). Also: a full sanitizer build can
regenerate `src/hb-ot-shaper-use-machine.hh` in the source tree; restore it,
never commit it.

Tools in `memory/gpu-tools/`: `n1_bound.py` (yield), `n1_validate.py`
(identity + bound validation), `zero_handles.py`, `hb_geo_single.py`.
