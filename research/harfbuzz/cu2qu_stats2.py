#!/usr/bin/env python3
# Temporary counters for the spline version of src/hb-gpu-cu2qu.hh (revert
# with git checkout).  Printed at exit when HB_CU2QU_STATS is set: curves
# per cubic, pieces-per-cubic histogram, trials, cubics no spline fitted,
# the largest equal-t deviation (16 samples per quadratic) and the tangent
# angles at joints and at the cubic's ends.
import sys
p = sys.argv[1]
s = open(p).read()

def rep(a, b):
    global s
    assert s.count(a) == 1, (s.count(a), a[:70])
    s = s.replace(a, b)

rep('''/*
 * A cubic in the power basis, to cut pieces out of.
 */
''', r'''#include <stdio.h>
struct hb_c2q_stats_t
{
  unsigned long long cubics = 0, quads = 0, hopeless = 0, trials = 0, first_try = 0;
  unsigned long long n_hist[130] = {0};
  double max_quad_dev = 0;
  double max_joint = 0, max_end = 0;
  unsigned long long joints = 0, joints_1 = 0, ends = 0, ends_1 = 0;
  void angle (double ax, double ay, double bx, double by, bool end)
  {
    double la = hypot (ax, ay), lb = hypot (bx, by);
    if (!(la > 1e-9) || !(lb > 1e-9)) return;
    double c = (ax * bx + ay * by) / (la * lb);
    if (c > 1) c = 1;
    if (c < -1) c = -1;
    double deg = acos (c) * 57.29577951308232;
    if (end) { ends++; if (deg > max_end) max_end = deg; if (deg > 1) ends_1++; }
    else { joints++; if (deg > max_joint) max_joint = deg; if (deg > 1) joints_1++; }
  }
  void quad (const double *c, double qx, double qy)
  {
    /* c: the piece's four control points; the quadratic runs from its
     * own end points, passed after them. */
    for (unsigned j = 1; j < 16; j++)
    {
      double t = j / 16., u = 1 - t;
      double bx = u*u*u*c[0] + 3*u*u*t*c[2] + 3*u*t*t*c[4] + t*t*t*c[6];
      double by = u*u*u*c[1] + 3*u*u*t*c[3] + 3*u*t*t*c[5] + t*t*t*c[7];
      double x = u*u*c[8] + 2*u*t*qx + t*t*c[10];
      double y = u*u*c[9] + 2*u*t*qy + t*t*c[11];
      double d = hypot (bx - x, by - y);
      if (d > max_quad_dev) max_quad_dev = d;
    }
  }
  ~hb_c2q_stats_t ()
  {
    if (!getenv ("HB_CU2QU_STATS")) return;
    fprintf (stderr, "cu2qu: cubics %llu quads %llu curves-per-cubic %.4f; no spline fitted %llu; trials per cubic %.3f; first trial fitted %.1f%%\n",
	     cubics, quads, cubics ? (double) quads / cubics : 0., hopeless,
	     cubics ? (double) trials / cubics : 0., cubics ? 100. * first_try / cubics : 0.);
    fprintf (stderr, "  max deviation: quadratics %.4f\n", max_quad_dev);
    fprintf (stderr, "  tangent angle: at %llu joints max %.3f deg (%llu over 1); at %llu cubic ends max %.3f deg (%llu over 1)\n",
	     joints, max_joint, joints_1, ends, max_end, ends_1);
    fprintf (stderr, "  pieces per cubic:");
    for (unsigned i = 0; i < 130; i++) if (n_hist[i]) fprintf (stderr, " %u:%llu", i, n_hist[i]);
    fprintf (stderr, "\n");
  }
};
static hb_c2q_stats_t hb_c2q_stats;
/*
 * A cubic in the power basis, to cut pieces out of.
 */
''')

rep('''  hb_gpu_cu2qu_point_t q1;
  if (hb_gpu_cubic_approx_quadratic (c0, c1, c2, c3, tolerance, &q1))
  {
    g->acc_conic_to (q1.x, q1.y, c3.x, c3.y);
    return;
  }
''', '''  hb_c2q_stats.cubics++;
  hb_gpu_cu2qu_point_t q1;
  if (hb_gpu_cubic_approx_quadratic (c0, c1, c2, c3, tolerance, &q1))
  {
    hb_c2q_stats.quads++; hb_c2q_stats.n_hist[1]++; hb_c2q_stats.trials++; hb_c2q_stats.first_try++;
    hb_c2q_stats.angle (q1.x - c0.x, q1.y - c0.y, c1.x - c0.x, c1.y - c0.y, true);
    hb_c2q_stats.angle (c3.x - q1.x, c3.y - q1.y, c3.x - c2.x, c3.y - c2.y, true);
    { double c[12] = {c0.x, c0.y, c1.x, c1.y, c2.x, c2.y, c3.x, c3.y, c0.x, c0.y, c3.x, c3.y};
      hb_c2q_stats.quad (c, q1.x, q1.y); }
    g->acc_conic_to (q1.x, q1.y, c3.x, c3.y);
    return;
  }
''')

rep('''  hb_gpu_cu2qu_point_t off[HB_GPU_CU2QU_MAX_N];
  hb_gpu_cu2qu_point_t on[HB_GPU_CU2QU_MAX_N + 1];
  for (; n <= HB_GPU_CU2QU_MAX_N; n++)
  {
    if (unlikely (!hb_budget_spend (*budget, HB_BUDGET_1, n)))
      return;

    if (!hb_gpu_cubic_approx_spline (cubic, n, tolerance, off, on))
      continue;
''', '''  hb_gpu_cu2qu_point_t off[HB_GPU_CU2QU_MAX_N];
  hb_gpu_cu2qu_point_t on[HB_GPU_CU2QU_MAX_N + 1];
  unsigned n_first = n;
  hb_c2q_stats.trials++;
  for (; n <= HB_GPU_CU2QU_MAX_N; n++)
  {
    if (unlikely (!hb_budget_spend (*budget, HB_BUDGET_1, n)))
      return;

    hb_c2q_stats.trials++;
    if (!hb_gpu_cubic_approx_spline (cubic, n, tolerance, off, on))
      continue;

    if (n == n_first) hb_c2q_stats.first_try++;
    hb_c2q_stats.quads += n; hb_c2q_stats.n_hist[n]++;
    hb_c2q_stats.angle (off[0].x - c0.x, off[0].y - c0.y, c1.x - c0.x, c1.y - c0.y, true);
    hb_c2q_stats.angle (c3.x - off[n - 1].x, c3.y - off[n - 1].y, c3.x - c2.x, c3.y - c2.y, true);
    for (unsigned k = 1; k < n; k++)
      hb_c2q_stats.angle (on[k].x - off[k - 1].x, on[k].y - off[k - 1].y, off[k].x - on[k].x, off[k].y - on[k].y, false);
    for (unsigned k = 0; k < n; k++)
    {
      hb_gpu_cu2qu_point_t pc[4];
      cubic.piece (k, n, pc);
      double c[12] = {pc[0].x, pc[0].y, pc[1].x, pc[1].y, pc[2].x, pc[2].y, pc[3].x, pc[3].y,
		      on[k].x, on[k].y, on[k + 1].x, on[k + 1].y};
      hb_c2q_stats.quad (c, off[k].x, off[k].y);
    }
''')

rep('''  g->acc_line_to (c1.x, c1.y);
  g->acc_line_to (c2.x, c2.y);
  g->acc_line_to (c3.x, c3.y);
}
''', '''  hb_c2q_stats.hopeless++;
  g->acc_line_to (c1.x, c1.y);
  g->acc_line_to (c2.x, c2.y);
  g->acc_line_to (c3.x, c3.y);
}
''')
open(p, 'w').write(s)
