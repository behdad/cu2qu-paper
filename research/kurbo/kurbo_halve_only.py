#!/usr/bin/env python3
"""Make the standalone "halve as cu2qu does" version of kurbo/src/cubicbez.rs:
upstream main plus only the halving fix (no closed-form n).

Usage: kurbo_halve_only.py <upstream cubicbez.rs> <final cubicbez.rs> <out>
"""
import sys

up = open(sys.argv[1]).read()
final = open(sys.argv[2]).read()


def block(text, start, end):
    i = text.index(start)
    j = text.index(end, i)
    return i, j


# 1. fit_inside (and subdivide_2 after it) from the final file
i, j = block(up, "    fn fit_inside(&self, distance: f64) -> bool {", "    /// Is this cubic Bezier curve finite?")
fi, fj = block(final, "    fn fit_inside(&self, distance: f64) -> bool {", "    /// Is this cubic Bezier curve finite?")
s = up[:i] + final[fi:fj] + up[j:]

# 2. split_into_n halves with subdivide_2
i, j = block(s, "    fn split_into_n(&self, n: usize) -> impl Iterator<Item = CubicBez> {", "    fn parameters(&self)")
body = s[i:j]
assert body.count(".subdivide()") == 5, body.count(".subdivide()")
s = s[:i] + body.replace(".subdivide()", ".subdivide_2()") + s[j:]

# 3. the end-point test of the spline, with the hypot that cu2qu uses
a = "            if d1.hypot() > accuracy\n"
assert s.count(a) == 1
s = s.replace(a, """            // `f64::hypot` and not `Vec2::hypot`, which rounds differently: at
            // the limit this has to decide as cu2qu does.
            if d1.x.hypot(d1.y) > accuracy
""")

# 4. tests, before the last test of the file
tests = '''    #[test]
    fn cubicbez_approx_spline_halves_as_cu2qu() {
        // Ensure a cubic is split in two with the arithmetic of fonttools cu2qu.
        // `ParamCurve::subdivide` rounds differently when coordinates are not
        // integers, and gave 1.199999999999999 for the x of the third point.
        let cubic = CubicBez::new((5.9, 3.6), (0.0, 9.2), (7.0, 6.6), (-16.2, -7.7));
        let spline = cubic.approx_spline(1.0).unwrap();
        assert_eq!(
            spline.points(),
            &[
                Point::new(5.9, 3.6),
                Point::new(1.4749999999999996, 7.799999999999999),
                Point::new(1.1999999999999993, 3.0250000000000012),
                Point::new(-16.2, -7.7),
            ]
        );
    }

    #[test]
    fn cubicbez_approx_spline_at_the_limit() {
        // The spline of six quadratics for this cubic is as far from it as the
        // accuracy, exactly. Testing it for fit used to split a piece without
        // end, and overflow the stack. cu2qu takes seven.
        let cubic = CubicBez::new(
            (386.0, -158.0),
            (527.0, -28.0),
            (480.0, -99.0),
            (2270.0, 2329.0),
        );
        assert_eq!(cubic.approx_spline(4.0).unwrap().points().len(), 7 + 2);

        // The single quadratic for these is as far from the cubic as the
        // accuracy, to the last bit, and its test did the same. cu2qu finds
        // that the single quadratic fits.
        for (cubic, accuracy) in [
            (
                CubicBez::new((-31.0, 97.0), (-51.0, 93.0), (-91.0, 76.0), (-139.0, 52.0)),
                1.987615979999813,
            ),
            (
                CubicBez::new((-20.0, 7.0), (-20.0, -27.5), (-12.0, -32.5), (4.0, -32.0)),
                3.5555555555555554,
            ),
        ] {
            assert_eq!(cubic.approx_spline(accuracy).unwrap().points().len(), 1 + 2);
        }
    }

    #[test]
    fn cubicbez_approx_spline_nan() {
        // These overflowed the stack too.
        let cubic = CubicBez::new((0.0, 0.0), (f64::NAN, 0.0), (1.0, 1.0), (2.0, 2.0));
        assert!(cubic.approx_spline(1.0).is_none());
        assert!(cubics_to_quadratic_splines(&[cubic], 1.0).is_none());

        let cubic = CubicBez::new((0.0, 0.0), (0.0, 100.0), (100.0, 200.0), (354.0, 300.0));
        assert!(cubic.approx_spline(f64::NAN).is_none());
        assert!(cubics_to_quadratic_splines(&[cubic], f64::NAN).is_none());
    }

'''
a = "    #[test]\n    fn cubics_to_quadratic_splines_matches_python() {"
assert s.count(a) == 1
s = s.replace(a, tests + a)
open(sys.argv[3], "w").write(s)
