#!/usr/bin/env python3
# Prototype V5 for src/hb-gpu-cu2qu.hh, applied to the pristine file (it
# applies V1 itself), with counters built in (HB_CU2QU_STATS; pass 'nostats'
# as a second argument to leave them out).
#
# The fonttools cu2qu spline, not its single-quadratic step applied
# recursively:
#   n = 1: the quadratic through the intersection of the end tangents, if it
#          is within tolerance (as now).
#   n > 1: split the cubic into n equal parameter intervals.  Piece k gets
#          the off-curve point
#            lerp (a + 1.5 (b - a), d + 1.5 (c - d), k / (n - 1))
#          of its own control points a b c d, and consecutive pieces meet at
#          the midpoint of their off-curve points.  The first and the last
#          off-curve points lie on the cubic's end tangents.  Each piece is
#          checked against its sub-cubic with the exact equal-t test, the
#          end-point offsets included.  Try n = 2, 3, ... until all fit.
import sys, subprocess, os
p = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
subprocess.check_call([sys.executable, os.path.join(here, 'c2q_fix_v1.py'), p])
s = open(p).read()
stats = len(sys.argv) < 3 or sys.argv[2] != 'nostats'

start = s.index('/*\n * Split a cubic at t = 0.5 (de Casteljau).')
end = s.index('#endif /* HB_GPU_CU2QU_HH */')
new = r'''#ifdef HB_C2Q_STATS
#include <stdio.h>
struct hb_c2q_stats_t
{
  unsigned long long cubics = 0, quads = 0, capped = 0, trials = 0;
  unsigned long long n_hist[66] = {0};
  double max_quad_dev = 0;
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
    fprintf (stderr, "cu2qu: cubics %llu quads %llu curves-per-cubic %.4f; cubics over the piece cap %llu; trials per cubic %.3f\n",
	     cubics, quads, cubics ? (double) quads / cubics : 0., capped, cubics ? (double) trials / cubics : 0.);
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
#define HB_GPU_CU2QU_MAX_PIECES 32

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
  if (unlikely (!g->success ||
		!hb_budget_spend (*g->get_budget_remaining (), HB_BUDGET_1)))
    return;

#ifdef HB_C2Q_STATS
  hb_c2q_stats.cubics++;
#endif

  hb_gpu_cu2qu_point_t q1;
  if (hb_gpu_cubic_approx_quadratic (c0, c1, c2, c3, tolerance, &q1))
  {
#ifdef HB_C2Q_STATS
    hb_c2q_stats.quads++;
    hb_c2q_stats.n_hist[1]++;
    hb_c2q_stats.trials++;
    hb_c2q_stats.angle (q1.x - c0.x, q1.y - c0.y, c1.x - c0.x, c1.y - c0.y, true);
    hb_c2q_stats.angle (c3.x - q1.x, c3.y - q1.y, c3.x - c2.x, c3.y - c2.y, true);
    for (unsigned k = 1; k < 16; k++)
    {
      double t = k / 16., u = 1 - t;
      double bx = u*u*u*c0.x + 3*u*u*t*c1.x + 3*u*t*t*c2.x + t*t*t*c3.x;
      double by = u*u*u*c0.y + 3*u*u*t*c1.y + 3*u*t*t*c2.y + t*t*t*c3.y;
      double qx = u*u*c0.x + 2*u*t*q1.x + t*t*c3.x;
      double qy = u*u*c0.y + 2*u*t*q1.y + t*t*c3.y;
      double d = hypot (bx - qx, by - qy);
      if (d > hb_c2q_stats.max_quad_dev) hb_c2q_stats.max_quad_dev = d;
    }
#endif
    g->acc_conic_to (q1.x, q1.y, c3.x, c3.y);
    return;
  }

  /* Power basis: B(t) = c0 + t (b1 + t (b2 + t b3)) */
  double b1x = 3 * (c1.x - c0.x), b1y = 3 * (c1.y - c0.y);
  double b2x = 3 * (c0.x - 2 * c1.x + c2.x), b2y = 3 * (c0.y - 2 * c1.y + c2.y);
  double b3x = c3.x - 3 * c2.x + 3 * c1.x - c0.x;
  double b3y = c3.y - 3 * c2.y + 3 * c1.y - c0.y;

  /* Off-curve point of each piece and the on-curve points between them. */
  hb_gpu_cu2qu_point_t off[HB_GPU_CU2QU_MAX_PIECES];
  hb_gpu_cu2qu_point_t on[HB_GPU_CU2QU_MAX_PIECES + 1];
  /* Sub-cubic control points, for the check. */
  hb_gpu_cu2qu_point_t sub[HB_GPU_CU2QU_MAX_PIECES][4];

  unsigned n = 2;
  for (;; n++)
  {
    if (unlikely (!hb_budget_spend (*g->get_budget_remaining (), HB_BUDGET_1, n)))
      return;
#ifdef HB_C2Q_STATS
    hb_c2q_stats.trials++;
#endif

    double dt = 1. / n;
    for (unsigned k = 0; k < n; k++)
    {
      double t0 = k * dt, t1 = k + 1 == n ? 1. : (k + 1) * dt;
      hb_gpu_cu2qu_point_t a = {c0.x + t0 * (b1x + t0 * (b2x + t0 * b3x)),
				c0.y + t0 * (b1y + t0 * (b2y + t0 * b3y))};
      hb_gpu_cu2qu_point_t d = {c0.x + t1 * (b1x + t1 * (b2x + t1 * b3x)),
				c0.y + t1 * (b1y + t1 * (b2y + t1 * b3y))};
      if (!k) a = c0;
      if (k + 1 == n) d = c3;
      double da_x = b1x + t0 * (2 * b2x + t0 * 3 * b3x), da_y = b1y + t0 * (2 * b2y + t0 * 3 * b3y);
      double dd_x = b1x + t1 * (2 * b2x + t1 * 3 * b3x), dd_y = b1y + t1 * (2 * b2y + t1 * 3 * b3y);
      hb_gpu_cu2qu_point_t b = {a.x + da_x * dt / 3, a.y + da_y * dt / 3};
      hb_gpu_cu2qu_point_t c = {d.x - dd_x * dt / 3, d.y - dd_y * dt / 3};
      sub[k][0] = a; sub[k][1] = b; sub[k][2] = c; sub[k][3] = d;

      hb_gpu_cu2qu_point_t p1 = {a.x + 1.5 * (b.x - a.x), a.y + 1.5 * (b.y - a.y)};
      hb_gpu_cu2qu_point_t p2 = {d.x + 1.5 * (c.x - d.x), d.y + 1.5 * (c.y - d.y)};
      off[k] = hb_gpu_cu2qu_lerp (p1, p2, (double) k / (n - 1));
    }
    on[0] = c0;
    for (unsigned k = 1; k < n; k++)
      on[k] = hb_gpu_cu2qu_lerp (off[k - 1], off[k], 0.5);
    on[n] = c3;

    if (n >= HB_GPU_CU2QU_MAX_PIECES)
    {
#ifdef HB_C2Q_STATS
      hb_c2q_stats.capped++;
#endif
      break;
    }

    bool fits = true;
    for (unsigned k = 0; k < n && fits; k++)
    {
      hb_gpu_cu2qu_point_t e0 = {on[k].x - sub[k][0].x, on[k].y - sub[k][0].y};
      hb_gpu_cu2qu_point_t e3 = {on[k + 1].x - sub[k][3].x, on[k + 1].y - sub[k][3].y};
      hb_gpu_cu2qu_point_t e1 = {on[k].x + (off[k].x - on[k].x) * (2.0 / 3.0) - sub[k][1].x,
				 on[k].y + (off[k].y - on[k].y) * (2.0 / 3.0) - sub[k][1].y};
      hb_gpu_cu2qu_point_t e2 = {on[k + 1].x + (off[k].x - on[k + 1].x) * (2.0 / 3.0) - sub[k][2].x,
				 on[k + 1].y + (off[k].y - on[k + 1].y) * (2.0 / 3.0) - sub[k][2].y};
      if (hypot (e3.x, e3.y) > tolerance ||
	  !hb_gpu_cubic_farthest_fit_inside (e0, e1, e2, e3, tolerance, 0))
	fits = false;
    }
    if (fits)
      break;
  }

#ifdef HB_C2Q_STATS
  hb_c2q_stats.n_hist[n]++;
  hb_c2q_stats.quads += n;
  hb_c2q_stats.angle (off[0].x - c0.x, off[0].y - c0.y, c1.x - c0.x, c1.y - c0.y, true);
  hb_c2q_stats.angle (c3.x - off[n - 1].x, c3.y - off[n - 1].y, c3.x - c2.x, c3.y - c2.y, true);
  for (unsigned k = 1; k < n; k++)
    hb_c2q_stats.angle (on[k].x - off[k - 1].x, on[k].y - off[k - 1].y, off[k].x - on[k].x, off[k].y - on[k].y, false);
  for (unsigned k = 0; k < n; k++)
    for (unsigned j = 1; j < 16; j++)
    {
      double t = j / 16., u = 1 - t;
      double bx = u*u*u*sub[k][0].x + 3*u*u*t*sub[k][1].x + 3*u*t*t*sub[k][2].x + t*t*t*sub[k][3].x;
      double by = u*u*u*sub[k][0].y + 3*u*u*t*sub[k][1].y + 3*u*t*t*sub[k][2].y + t*t*t*sub[k][3].y;
      double qx = u*u*on[k].x + 2*u*t*off[k].x + t*t*on[k + 1].x;
      double qy = u*u*on[k].y + 2*u*t*off[k].y + t*t*on[k + 1].y;
      double d = hypot (bx - qx, by - qy);
      if (d > hb_c2q_stats.max_quad_dev) hb_c2q_stats.max_quad_dev = d;
    }
#endif

  for (unsigned k = 0; k < n; k++)
  {
    if (unlikely (hb_gpu_cu2qu_unencodable (on[k + 1])))
    {
      g->acc_line_to (on[k + 1].x, on[k + 1].y);
      g->acc_line_to (c3.x, c3.y);
      return;
    }
    g->acc_conic_to (off[k].x, off[k].y, on[k + 1].x, on[k + 1].y);
  }
}


'''
if stats:
    new = '#define HB_C2Q_STATS 1\n' + new
s = s[:start] + new + s[end:]
open(p, 'w').write(s)
