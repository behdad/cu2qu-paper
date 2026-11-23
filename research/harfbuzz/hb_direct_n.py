#!/usr/bin/env python3
# Applies the "n in closed form" change to src/hb-gpu-cu2qu.hh (main as of
# 0c72f9596): the spline's error is known from n and the third difference,
# so n is computed, and the trial loop and the per-piece test go.
import sys
p = sys.argv[1]
s = open(p).read()

def rep(a, b):
    global s
    assert s.count(a) == 1, (s.count(a), a[:70])
    s = s.replace(a, b)

start = s.index('/*\n * Try to approximate a cubic Bézier with a spline of n > 1 quadratics.')
end = s.index('#endif /* HB_GPU_CU2QU_HH */')
new = '''/*
 * The spline of n > 1 quadratics for a cubic.
 *
 * The cubic is cut into n equal pieces.  Each gives one off-curve point,
 * between the points its two end tangents suggest, nearer the first one
 * in the first pieces and the second one in the last, and consecutive
 * quadratics meet halfway between their off-curve points.  So the spline
 * is smooth, and it leaves both ends of the cubic along the cubic's own
 * tangents.
 */

/* The off-curve point for the k-th of n pieces, @p its control points. */
static inline hb_gpu_cu2qu_point_t
hb_gpu_cubic_spline_control (const hb_gpu_cu2qu_point_t p[4],
			     unsigned k, unsigned n)
{
  hb_gpu_cu2qu_point_t q0 = {p[0].x + (p[1].x - p[0].x) * 1.5,
			     p[0].y + (p[1].y - p[0].y) * 1.5};
  hb_gpu_cu2qu_point_t q1 = {p[3].x + (p[2].x - p[3].x) * 1.5,
			     p[3].y + (p[2].y - p[3].y) * 1.5};
  return hb_gpu_cu2qu_lerp (q0, q1, (double) k / (n - 1));
}

/*
 * How far the spline of n > 1 quadratics gets from its cubic at equal
 * parameters, as a multiple of the cubic's third difference over n^3.
 *
 * It does not depend on the cubic.  With d the third difference over
 * n^3, the off-curve point of piece k is the one of least error for the
 * piece plus d (1 - 2k/(n-1)) / 4, which puts the on-curve point between
 * pieces k-1 and k at d g(k) from the cubic, g(k) = (n-2k) / (4(n-1)).
 * The difference between quadratic k and its piece then is d times
 *
 *   g(k) - u (2g(k) + k/(n-1)) + u^2 (1 + g(k) + g(k+1) + k/(n-1)) - u^3
 *
 * at parameter u, with g(0) = g(n) = 0.  It is largest in the first
 * piece or in the second, and likewise in the last two.
 */
static inline double
hb_gpu_cubic_spline_error (unsigned n)
{
  if (n == 2)
    return 4. / 27.;

  double g = (n - 2) / (4. * (n - 1));
  double first = 4. / 27. * (1 + g) * (1 + g) * (1 + g);

  double w2 = (.25 - .5 / (n - 1)) / 3;
  double second = (n - 3) / (4. * (n - 1)) + 2 * w2 * sqrt (w2);

  return hb_max (first, second);
}


/*
 * Convert a cubic into quadratic segments.
 */
static void
hb_gpu_cubic_to_quadratics (hb_gpu_draw_t *g,
			    hb_gpu_cu2qu_point_t c0,
			    hb_gpu_cu2qu_point_t c1,
			    hb_gpu_cu2qu_point_t c2,
			    hb_gpu_cu2qu_point_t c3,
			    double tolerance)
{
  /* Charge the work budget a unit for the fit test and one for each
   * quadratic computed, and stop once the accumulator has failed or the
   * budget is spent: whatever we would emit from here is dropped by
   * acc_emit() anyway. */
  int64_t *budget = g->get_budget_remaining ();
  if (unlikely (!g->success ||
		!hb_budget_spend (*budget, HB_BUDGET_1)))
    return;

  hb_gpu_cu2qu_point_t q1;
  if (hb_gpu_cubic_approx_quadratic (c0, c1, c2, c3, tolerance, &q1))
  {
    g->acc_conic_to (q1.x, q1.y, c3.x, c3.y);
    return;
  }

  hb_gpu_cu2qu_cubic_t cubic (c0, c1, c2, c3);

  /* The least n whose spline is within tolerance.  Written so that NaN
   * goes past the last. */
  double d3 = hypot (cubic.b3.x, cubic.b3.y);
  unsigned n = 2;
  while (n <= HB_GPU_CU2QU_MAX_N &&
	 !(d3 * hb_gpu_cubic_spline_error (n) <= tolerance * n * n * n))
    n++;

  if (unlikely (n > HB_GPU_CU2QU_MAX_N))
  {
    /* Too large for that many, which is far beyond what the encoder can
     * store, as when a transform puts an outline astronomically far
     * away: hb_gpu_draw_encode() will reject the outline.  Emit the
     * control polygon, which has the cubic's points in it. */
    g->acc_line_to (c1.x, c1.y);
    g->acc_line_to (c2.x, c2.y);
    g->acc_line_to (c3.x, c3.y);
    return;
  }

  if (unlikely (!hb_budget_spend (*budget, HB_BUDGET_1, n)))
    return;

  hb_gpu_cu2qu_point_t piece[4];
  cubic.piece (0, n, piece);
  hb_gpu_cu2qu_point_t off = hb_gpu_cubic_spline_control (piece, 0, n);
  for (unsigned k = 1; k < n; k++)
  {
    cubic.piece (k, n, piece);
    hb_gpu_cu2qu_point_t next_off = hb_gpu_cubic_spline_control (piece, k, n);
    hb_gpu_cu2qu_point_t on = hb_gpu_cu2qu_lerp (off, next_off, 0.5);
    g->acc_conic_to (off.x, off.y, on.x, on.y);
    off = next_off;
  }
  g->acc_conic_to (off.x, off.y, c3.x, c3.y);
}


'''
s = s[:start] + new + s[end:]

rep('''/* Most quadratics one cubic is converted to.  Glyphs take up to 10 at
 * this tolerance, and the wildest cubics whose control points the encoder
 * can store about 35.  One that would take more is emitted as lines. */''',
'''/* Most quadratics one cubic is converted to.  Glyphs take up to 10 at
 * this tolerance, and the wildest cubics whose control points the encoder
 * can store about 35.  One that would take more is emitted as lines. */''')
open(p, 'w').write(s)
