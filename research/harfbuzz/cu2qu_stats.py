#!/usr/bin/env python3
# Temporary instrumentation of src/hb-gpu-cu2qu.hh (reverted with git
# checkout).  Printed at exit when HB_CU2QU_STATS is set:
#  - cubics, nodes, quadratics and cap lines emitted;
#  - why the single-quadratic fit failed at each node and at the pieces that
#    reach the depth cap;
#  - the largest distance between an emitted quadratic and its cubic piece
#    at equal parameters (16 samples), and between a cap line and its piece.
# Works on the file before and after the fix: it only needs the anchors.
import sys
p = sys.argv[1]
s = open(p).read()

def rep(a, b, count=1):
    global s
    assert s.count(a) >= 1, a
    s = s.replace(a, b, count)

rep('''/*
 * Try to approximate a cubic Bézier with a single quadratic.
 */
''', '''#include <stdio.h>
struct hb_c2q_stats_t
{
  unsigned long long cubics = 0, nodes = 0, quads = 0, cap_lines = 0;
  /* reasons: 1 start tangent zero, 2 end tangent zero, 3 both zero,
   * 4 tangents parallel, 5 error over tolerance */
  unsigned long long fail[6] = {0}, cap[6] = {0};
  unsigned long long depth_hist[12] = {0};
  double max_quad_dev = 0, max_line_dev = 0;
  int last_reason = 0;
  ~hb_c2q_stats_t ()
  {
    if (!getenv ("HB_CU2QU_STATS")) return;
    fprintf (stderr, "cu2qu: cubics %llu nodes %llu quads %llu cap-lines %llu curves-per-cubic %.4f\\n",
	     cubics, nodes, quads, cap_lines, cubics ? (double) (quads + cap_lines) / cubics : 0.);
    fprintf (stderr, "  fit failed at nodes: start-zero %llu end-zero %llu both-zero %llu parallel %llu error %llu\\n",
	     fail[1], fail[2], fail[3], fail[4], fail[5]);
    fprintf (stderr, "  reached the cap:     start-zero %llu end-zero %llu both-zero %llu parallel %llu error %llu\\n",
	     cap[1], cap[2], cap[3], cap[4], cap[5]);
    fprintf (stderr, "  max deviation: quadratics %.4f, cap lines %.4f\\n", max_quad_dev, max_line_dev);
    fprintf (stderr, "  leaf depth:");
    for (unsigned i = 0; i < 12; i++) if (depth_hist[i]) fprintf (stderr, " %u:%llu", i, depth_hist[i]);
    fprintf (stderr, "\\n");
  }
};
static hb_c2q_stats_t hb_c2q_stats;
/*
 * Try to approximate a cubic Bézier with a single quadratic.
 */
''')

# Record why the fit failed.  The classification is computed on entry from
# the control points, so it is the same before and after the fix.
rep('''			       double tolerance,
			       hb_gpu_cu2qu_point_t *q1)
{
''', '''			       double tolerance,
			       hb_gpu_cu2qu_point_t *q1)
{
  {
    bool sz = c1.x == c0.x && c1.y == c0.y;
    bool ez = c3.x == c2.x && c3.y == c2.y;
    double cr = -(c1.y - c0.y) * (c3.x - c2.x) + (c1.x - c0.x) * (c3.y - c2.y);
    hb_c2q_stats.last_reason = sz && ez ? 3 : sz ? 1 : ez ? 2 : fabs (cr) < 1e-12 ? 4 : 5;
  }
''')

rep('''  hb_gpu_cu2qu_point_t q1;
  if (hb_gpu_cubic_approx_quadratic (c0, c1, c2, c3, tolerance, &q1))
  {
    g->acc_conic_to (q1.x, q1.y, c3.x, c3.y);
    return;
  }

  if (depth >= HB_GPU_CU2QU_MAX_DEPTH)
  {
    g->acc_line_to (c3.x, c3.y);
    return;
  }
''', '''  if (!depth) hb_c2q_stats.cubics++;
  hb_c2q_stats.nodes++;
  hb_gpu_cu2qu_point_t q1;
  if (hb_gpu_cubic_approx_quadratic (c0, c1, c2, c3, tolerance, &q1))
  {
    hb_c2q_stats.quads++;
    hb_c2q_stats.depth_hist[depth]++;
    for (unsigned i = 1; i < 16; i++)
    {
      double t = i / 16., u = 1 - t;
      double bx = u*u*u*c0.x + 3*u*u*t*c1.x + 3*u*t*t*c2.x + t*t*t*c3.x;
      double by = u*u*u*c0.y + 3*u*u*t*c1.y + 3*u*t*t*c2.y + t*t*t*c3.y;
      double qx = u*u*c0.x + 2*u*t*q1.x + t*t*c3.x;
      double qy = u*u*c0.y + 2*u*t*q1.y + t*t*c3.y;
      double d = hypot (bx - qx, by - qy);
      if (d > hb_c2q_stats.max_quad_dev) hb_c2q_stats.max_quad_dev = d;
    }
    g->acc_conic_to (q1.x, q1.y, c3.x, c3.y);
    return;
  }
  hb_c2q_stats.fail[hb_c2q_stats.last_reason]++;

  if (depth >= HB_GPU_CU2QU_MAX_DEPTH)
  {
    hb_c2q_stats.cap_lines++;
    hb_c2q_stats.cap[hb_c2q_stats.last_reason]++;
    hb_c2q_stats.depth_hist[depth]++;
    for (unsigned i = 1; i < 16; i++)
    {
      double t = i / 16., u = 1 - t;
      double bx = u*u*u*c0.x + 3*u*u*t*c1.x + 3*u*t*t*c2.x + t*t*t*c3.x;
      double by = u*u*u*c0.y + 3*u*u*t*c1.y + 3*u*t*t*c2.y + t*t*t*c3.y;
      double lx = u*c0.x + t*c3.x, ly = u*c0.y + t*c3.y;
      double d = hypot (bx - lx, by - ly);
      if (d > hb_c2q_stats.max_line_dev) hb_c2q_stats.max_line_dev = d;
    }
    g->acc_line_to (c3.x, c3.y);
    return;
  }
''')
open(p, 'w').write(s)
