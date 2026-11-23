#!/usr/bin/env python3
"""Stack the "compute n" change on the halving commit of kurbo.

Usage: kurbo_compute_n_on_halve.py <halve cubicbez.rs> <earlier final cubicbez.rs> <out>

Takes the pieces of the computed-n change from the earlier combined version
and applies them to the halving commit's file. Distances use kurbo's own
Vec2::hypot throughout.
"""
import sys

s = open(sys.argv[1]).read()
final = open(sys.argv[2]).read()


def between(text, start, end):
    i = text.index(start)
    j = text.index(end, i)
    return text[i:j]


def replace_between(text, start, end, new):
    i = text.index(start)
    j = text.index(end, i)
    return text[:i] + new + text[j:]


def rep(text, a, b):
    assert text.count(a) == 1, (text.count(a), a[:70])
    return text.replace(a, b)


# 1. the constant
const = between(final, "/// A single quadratic is not tried for a third difference", "/// A single cubic Bézier segment.")
s = rep(s, "const MAX_SPLINE_SPLIT: usize = 100;\n\n", "const MAX_SPLINE_SPLIT: usize = 100;\n\n" + const)

# 2. approx_spline, approx_spline_n, third_difference
blk = between(final, "    pub fn approx_spline(&self, accuracy: f64) -> Option<QuadSpline> {", "    fn approx_quad_control(&self, t: f64) -> Point {")
i = blk.index("    /// The length of the third difference of the control points.")
blk = blk[:i] + '''    /// The length of the third difference of the control points.
    ///
    /// Zero for a cubic that is a quadratic, and what decides how many
    /// quadratics another one takes.
    fn third_difference(&self) -> f64 {
        (self.p3.to_vec2() - 3.0 * self.p2.to_vec2() + 3.0 * self.p1.to_vec2() - self.p0.to_vec2())
            .hypot()
    }

'''
s = replace_between(s, "    pub fn approx_spline(&self, accuracy: f64) -> Option<QuadSpline> {", "    fn approx_quad_control(&self, t: f64) -> Point {", blk)

# 3. the single quadratic is not tried when the third difference rules it out
s = rep(s, """    fn try_approx_quadratic(&self, accuracy: f64) -> Option<QuadBez> {
        if let Some(q1) = Line::new(self.p0, self.p1)""", """    fn try_approx_quadratic(&self, accuracy: f64) -> Option<QuadBez> {
        // Not worth trying if the third difference rules it out. (A NaN does too.)
        let may_fit = self.third_difference() <= accuracy * QUADRATIC_MAX_D3;
        if !may_fit {
            return None;
        }
        if let Some(q1) = Line::new(self.p0, self.p1)""")

# 4. cubics_to_quadratic_splines, spline_error_factor, spline_count
blk = between(final, "pub fn cubics_to_quadratic_splines(", "#[cfg(test)]")
blk = rep(blk, """/// Rust port of cu2qu `spline_error_factor`, with the same arithmetic for the
/// same result to the last bit.
""", """/// Rust port of cu2qu `spline_error_factor`.
""")
s = replace_between(s, "pub fn cubics_to_quadratic_splines(", "#[cfg(test)]", blk)

# 5. tests
s = rep(s, """mod tests {
    use crate::{""", """mod tests {
    use super::{MAX_SPLINE_SPLIT, spline_error_factor};
    use crate::{""")
s = rep(s, """        let spline = c1.approx_spline_n(2, 343.0);
        assert!(spline.is_some());""", """        let spline = c1.approx_spline(343.0);
        assert!(spline.is_some());""")
s = rep(s, between(s, "        // FontTools can approximate this curve successfully in 7 splits, we can too", "    #[test]\n    fn cubic_to_quadratic_all_points_equal() {"),
        between(final, "        // FontTools approximates this curve with 6 quadratics, we do too", "    /// The largest distance between a cubic and its spline at equal parameters."))
new_tests = between(final, "    /// The largest distance between a cubic and its spline at equal parameters.", "    #[test]\n    fn cubicbez_approx_spline_third_difference_at_the_limit() {")
new_tests += '''    #[test]
    fn cubicbez_approx_spline_third_difference_at_the_limit() {
        // Two quadratics take a third difference of up to 54 times the accuracy,
        // three of up to 128 times: exactly that much still fits.
        for (x, n) in [(354.0, 2), (355.0, 3), (428.0, 3), (429.0, 4)] {
            let cubic = CubicBez::new((0.0, 0.0), (0.0, 100.0), (100.0, 200.0), (x, 300.0));
            assert_eq!(cubic.approx_spline(1.0).unwrap().points().len(), n + 2);
        }
    }

'''
new_tests += between(final, "    #[test]\n    fn cubicbez_approx_spline_more_quadratics_fit_too() {", "    #[test]\n    fn cubic_to_quadratic_all_points_equal() {")
new_tests = rep(new_tests, "        assert_eq!(d3, 270.0_f64.hypot(10.0));\n", "        assert_eq!(d3, 73000.0_f64.sqrt());\n")
s = rep(s, "    #[test]\n    fn cubic_to_quadratic_all_points_equal() {", new_tests + "    #[test]\n    fn cubic_to_quadratic_all_points_equal() {")
assert "x.hypot(" not in s, "libm hypot left in"
open(sys.argv[3], "w").write(s)
