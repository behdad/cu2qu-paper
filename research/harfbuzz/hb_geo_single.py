#!/usr/bin/env python3
"""Replace hb-gpu's recursive single-quadratic test with the closed-form
geometric bound.

Usage: hb_geo_single.py <path to src/hb-gpu-cu2qu.hh>
"""
import sys

p = sys.argv[1]
s = open(p).read()

i = s.index("/*\n * Check whether a cubic error curve stays within `tolerance` of the origin.\n */")
j = s.index("/*\n * A cubic in the power basis, to cut pieces out of.\n */")
new = r'''/*
 * Try to approximate a cubic Bézier with a single quadratic: the one
 * whose control point q is where the cubic's two end tangents meet.
 *
 * Returns whether the two curves are within `tolerance` of each other as
 * sets of points, not at equal parameters: the cubic runs along that
 * quadratic at another pace, and what counts is how far it leaves it.
 *
 * With a = q - c0 and z = c3 - q, the cubic's control points are at
 * c1 = c0 + l1 a and c2 = c3 - l2 z, and the cubic B is the quadratic Q
 * at another parameter plus an offset along the parabola's axis:
 *
 *   B(t) = Q(tau(t)) + rho(t) (a - z)
 *
 *   tau(t) = 3 t (1-t)^2 l1/2 + 3 t^2 (1-t) (1 - l2/2) + t^3
 *   rho(t) = 3 t^2 (1-t)^2 (l1 + l2 - 4/3 + 3 e(t)^2)
 *   e(t)   = (1-t) (l1/2 - 1/3) + t (1/3 - l2/2)
 *
 * which is an identity of polynomials.  For l1 and l2 from 0 to 2, tau
 * is a Bézier function with control values from 0 to 1, so it takes
 * [0,1] onto [0,1]: every point of the cubic is |rho| |a - z| from a
 * point of the quadratic, and every point of the quadratic that far from
 * one of the cubic.  As t^2 (1-t)^2 is at most 1/16, the two are within
 *
 *   3/16 max |l1 + l2 - 4/3 + 3 e^2| |a - z|
 *
 * So a single quadratic fits when the two handles add up to 4/3 of the
 * way to where the tangents meet, however that is split between them.
 * (At equal parameters it would take 2/3 each.)
 */
static bool
hb_gpu_cubic_approx_quadratic (hb_gpu_cu2qu_point_t c0,
			       hb_gpu_cu2qu_point_t c1,
			       hb_gpu_cu2qu_point_t c2,
			       hb_gpu_cu2qu_point_t c3,
			       double tolerance,
			       hb_gpu_cu2qu_point_t *q1)
{
  double ax = c1.x - c0.x, ay = c1.y - c0.y;
  double dx = c3.x - c2.x, dy = c3.y - c2.y;

  /* A control point on its end point leaves that end no tangent to
   * intersect.  The curve leaves the end toward the other control point,
   * so take the direction from that one: the quadratic's control point
   * is then that control point. */
  if (ax == 0 && ay == 0)
  {
    ax = c2.x - c0.x;
    ay = c2.y - c0.y;
  }
  if (dx == 0 && dy == 0)
  {
    dx = c3.x - c1.x;
    dy = c3.y - c1.y;
  }

  double px = -ay, py = ax;

  double denom = px * dx + py * dy;
  if (fabs (denom) < 1e-12)
    return false;

  double h = (px * (c0.x - c2.x) + py * (c0.y - c2.y)) / denom;

  q1->x = c2.x + dx * h;
  q1->y = c2.y + dy * h;

  /* The two legs of the quadratic, and how far along them the cubic's
   * control points are.  Written so that NaN gives up. */
  hb_gpu_cu2qu_point_t a = {q1->x - c0.x, q1->y - c0.y};
  hb_gpu_cu2qu_point_t z = {c3.x - q1->x, c3.y - q1->y};
  double aa = a.x * a.x + a.y * a.y;
  double zz = z.x * z.x + z.y * z.y;
  if (!(aa > 0 && zz > 0))
    return false;
  double l1 = ((c1.x - c0.x) * a.x + (c1.y - c0.y) * a.y) / aa;
  double l2 = ((c3.x - c2.x) * z.x + (c3.y - c2.y) * z.y) / zz;
  if (!(l1 >= 0 && l1 <= 2 && l2 >= 0 && l2 <= 2))
    return false;

  /* e(t) goes from p to r, so e^2 is between the smaller of their
   * squares, or zero if they differ in sign, and the larger. */
  double p = l1 * .5 - 1. / 3., r = 1. / 3. - l2 * .5;
  double e2_max = hb_max (p * p, r * r);
  double e2_min = p * r <= 0 ? 0. : hb_min (p * p, r * r);
  double g = l1 + l2 - 4. / 3.;
  double m = hb_max (fabs (g + 3 * e2_max), fabs (g + 3 * e2_min));

  return 3. / 16. * m * hypot (a.x - z.x, a.y - z.y) <= tolerance;
}


'''
s = s[:i] + new + s[j:]

a = """ * Cubic-to-quadratic Bézier conversion.
 * Ported from fonttools cu2qu:
 * https://github.com/fonttools/fonttools/blob/main/Lib/fontTools/cu2qu/cu2qu.py
"""
assert s.count(a) == 1
open(p, "w").write(s)
