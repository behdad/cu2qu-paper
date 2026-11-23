#!/usr/bin/env python3
# Prototype V3 for src/hb-gpu-cu2qu.hh, applied to the pristine file, with
# the same counters as cu2qu_stats.py built in (HB_CU2QU_STATS).
#
# Closed-form conversion: split the cubic into n equal parameter intervals
# and replace each piece by the quadratic whose control point is
# (3 (p1 + p2) - p0 - p3) / 4.  The difference between a cubic and that
# quadratic is (d3 / 2) t (1 - t) (2 t - 1), d3 the cubic's third
# difference, so it is at most sqrt(3)/36 |d3|, and a piece of 1/n of the
# cubic has d3 / n^3: n = ceil (cbrt (sqrt(3)/36 |d3| / tolerance)).
import sys
p = sys.argv[1]
s = open(p).read()
start = s.index('/*\n * Check whether a cubic error curve stays within')
end = s.index('#endif /* HB_GPU_CU2QU_HH */')
new = r'''#ifdef HB_C2Q_STATS
#include <stdio.h>
struct hb_c2q_stats_t
{
  unsigned long long cubics = 0, nodes = 0, quads = 0, capped = 0;
  unsigned long long n_hist[66] = {0};
  double max_quad_dev = 0;
  /* Angles in degrees: between the two quadratics at a joint within a
   * cubic, and between a quadratic and the cubic's own tangent at the
   * cubic's ends. */
  double max_joint = 0, max_end = 0;
  unsigned long long joints = 0, joints_1 = 0, joints_5 = 0, ends = 0, ends_1 = 0, ends_5 = 0;
  double sum_joint = 0, sum_end = 0;
  void angle (double ax, double ay, double bx, double by, bool end)
  {
    double la = hypot (ax, ay), lb = hypot (bx, by);
    if (!(la > 1e-9) || !(lb > 1e-9)) return;
    double c = (ax * bx + ay * by) / (la * lb);
    if (c > 1) c = 1;
    if (c < -1) c = -1;
    double deg = acos (c) * 57.29577951308232;
    if (end)
    { ends++; sum_end += deg; if (deg > max_end) max_end = deg; if (deg > 1) ends_1++; if (deg > 5) ends_5++; }
    else
    { joints++; sum_joint += deg; if (deg > max_joint) max_joint = deg; if (deg > 1) joints_1++; if (deg > 5) joints_5++; }
  }
  ~hb_c2q_stats_t ()
  {
    if (!getenv ("HB_CU2QU_STATS")) return;
    fprintf (stderr, "cu2qu: cubics %llu quads %llu curves-per-cubic %.4f; cubics over the piece cap %llu\n",
	     cubics, quads, cubics ? (double) quads / cubics : 0., capped);
    fprintf (stderr, "  max deviation: quadratics %.4f\n", max_quad_dev);
    fprintf (stderr, "  kinks at joints within a cubic: %llu, mean %.3f deg, max %.2f, over 1 deg %llu, over 5 deg %llu\n",
	     joints, joints ? sum_joint / joints : 0., max_joint, joints_1, joints_5);
    fprintf (stderr, "  tangent off at cubic ends:      %llu, mean %.3f deg, max %.2f, over 1 deg %llu, over 5 deg %llu\n",
	     ends, ends ? sum_end / ends : 0., max_end, ends_1, ends_5);
    fprintf (stderr, "  pieces per cubic:");
    for (unsigned i = 0; i < 66; i++) if (n_hist[i]) fprintf (stderr, " %u:%llu", i, n_hist[i]);
    fprintf (stderr, "\n");
  }
};
static hb_c2q_stats_t hb_c2q_stats;
#endif

/* Most pieces one cubic is split into. */
#define HB_GPU_CU2QU_MAX_PIECES 64

/*
 * Convert a cubic into quadratic segments.
 */
static void
hb_gpu_cubic_to_quadratics (hb_gpu_draw_t *g,
			    hb_gpu_cu2qu_point_t c0,
			    hb_gpu_cu2qu_point_t c1,
			    hb_gpu_cu2qu_point_t c2,
			    hb_gpu_cu2qu_point_t c3,
			    double tolerance, unsigned depth HB_UNUSED)
{
  /* Third difference. */
  double ex = c3.x - 3 * c2.x + 3 * c1.x - c0.x;
  double ey = c3.y - 3 * c2.y + 3 * c1.y - c0.y;
  double err2 = ex * ex + ey * ey;

  /* n^6 >= err2 / (432 tolerance^2) */
  double scale = 1. / (432. * tolerance * tolerance);
  double want = err2 * scale;
  unsigned n = 1;
  while (n < HB_GPU_CU2QU_MAX_PIECES && !((double) n * n * n * n * n * n >= want))
    n++;
#ifdef HB_C2Q_STATS
  hb_c2q_stats.cubics++;
  hb_c2q_stats.n_hist[n]++;
  if (!((double) n * n * n * n * n * n >= want)) hb_c2q_stats.capped++;
  double prev_qx = 0, prev_qy = 0;
#endif

  /* Power basis: B(t) = c0 + t (b1 + t (b2 + t b3)) */
  double b1x = 3 * (c1.x - c0.x), b1y = 3 * (c1.y - c0.y);
  double b2x = 3 * (c0.x - 2 * c1.x + c2.x), b2y = 3 * (c0.y - 2 * c1.y + c2.y);
  double b3x = ex, b3y = ey;

  double dt = 1. / n;
  double p0x = c0.x, p0y = c0.y;
  double d0x = b1x, d0y = b1y; /* B'(0) */
  for (unsigned i = 1; i <= n; i++)
  {
    if (unlikely (!g->success ||
		  !hb_budget_spend (*g->get_budget_remaining (), HB_BUDGET_1)))
      return;

    double t = i == n ? 1. : i * dt;
    double p3x, p3y;
    if (i == n) { p3x = c3.x; p3y = c3.y; }
    else
    {
      p3x = c0.x + t * (b1x + t * (b2x + t * b3x));
      p3y = c0.y + t * (b1y + t * (b2y + t * b3y));
    }
    double d1x = b1x + t * (2 * b2x + t * 3 * b3x);
    double d1y = b1y + t * (2 * b2y + t * 3 * b3y);

    /* (3 (p1 + p2) - p0 - p3) / 4 with p1 = p0 + B'(t0) dt / 3 and
     * p2 = p3 - B'(t1) dt / 3. */
    double qx = (p0x + p3x) * .5 + (d0x - d1x) * (dt * .25);
    double qy = (p0y + p3y) * .5 + (d0y - d1y) * (dt * .25);

#ifdef HB_C2Q_STATS
    hb_c2q_stats.quads++;
    if (i == 1)
      hb_c2q_stats.angle (qx - p0x, qy - p0y, c1.x - c0.x, c1.y - c0.y, true);
    else
      hb_c2q_stats.angle (p0x - prev_qx, p0y - prev_qy, qx - p0x, qy - p0y, false);
    if (i == n)
      hb_c2q_stats.angle (p3x - qx, p3y - qy, c3.x - c2.x, c3.y - c2.y, true);
    prev_qx = qx; prev_qy = qy;
    {
      /* Deviation from the cubic piece at equal parameters. */
      double t0 = (i - 1) * dt;
      for (unsigned k = 1; k < 16; k++)
      {
	double u = k / 16., tt = t0 + u * dt;
	double bx = c0.x + tt * (b1x + tt * (b2x + tt * b3x));
	double by = c0.y + tt * (b1y + tt * (b2y + tt * b3y));
	double qqx = (1-u)*(1-u)*p0x + 2*u*(1-u)*qx + u*u*p3x;
	double qqy = (1-u)*(1-u)*p0y + 2*u*(1-u)*qy + u*u*p3y;
	double d = hypot (bx - qqx, by - qqy);
	if (d > hb_c2q_stats.max_quad_dev) hb_c2q_stats.max_quad_dev = d;
      }
    }
#endif

    if (unlikely (hb_gpu_cu2qu_unencodable ({p3x, p3y})))
    {
      /* As before: the outline cannot be encoded; finish with lines. */
      g->acc_line_to (p3x, p3y);
      g->acc_line_to (c3.x, c3.y);
      return;
    }

    g->acc_conic_to (qx, qy, p3x, p3y);

    p0x = p3x; p0y = p3y;
    d0x = d1x; d0y = d1y;
  }
}


'''
if len(sys.argv) < 3 or sys.argv[2] != 'nostats':
    new = '#define HB_C2Q_STATS 1\n' + new
s = s[:start] + new + s[end:]
open(p, 'w').write(s)
