#!/usr/bin/env python3
"""From the final kurbo/src/cubicbez.rs, make the state of the FIRST commit
(compute the number of quadratics) by backing out what belongs to the second
(halve as cu2qu does; libm hypot in the fit test).

Usage: kurbo_split_commits.py <final cubicbez.rs> <out>
"""
import re, sys

s = open(sys.argv[1]).read()


def cut(text, start, end):
    """Remove from the line containing `start` up to, not including, `end`."""
    i = text.index(start)
    j = text.index(end, i)
    return text[:i] + text[j:]


# 1. fit_inside back to upstream, subdivide_2 gone
i = s.index("    fn fit_inside(&self, distance: f64) -> bool {")
j = s.index("    /// Is this cubic Bezier curve finite?")
upstream_fit_inside = '''    fn fit_inside(&self, distance: f64) -> bool {
        if self.p2.to_vec2().hypot() <= distance && self.p1.to_vec2().hypot() <= distance {
            return true;
        }
        let mid =
            (self.p0.to_vec2() + 3.0 * (self.p1.to_vec2() + self.p2.to_vec2()) + self.p3.to_vec2())
                * 0.125;
        if mid.hypot() > distance {
            return false;
        }
        // Split in two. Note that cu2qu here uses a 3/8 subdivision. I don't know why.
        let (left, right) = self.subdivide();
        left.fit_inside(distance) && right.fit_inside(distance)
    }

'''
s = s[:i] + upstream_fit_inside + s[j:]

# 2. split_into_n back to subdivide()
assert s.count(".subdivide_2()") == 5, s.count(".subdivide_2()")
s = s.replace(".subdivide_2()", ".subdivide()")

# 3. the two tests of the second commit
s = cut(s, "    #[test]\n    fn cubicbez_approx_spline_single_quadratic_at_the_limit() {", "    #[test]\n    fn cubicbez_approx_spline_halves_as_cu2qu() {")
i = s.index("    #[test]\n    fn cubicbez_approx_spline_halves_as_cu2qu() {")
j = s.index("    #[test]", i + 10)
s = s[:i] + s[j:]
assert "subdivide_2" not in s and "halves_as_cu2qu" not in s and "single_quadratic_at_the_limit" not in s
open(sys.argv[2], "w").write(s)
