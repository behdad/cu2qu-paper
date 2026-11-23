/* gt.c -- measuring instrument ("ground truth") for the distance between a
 * cubic Bezier and the quadratic spline cu2qu-style conversion gives for it.
 *
 * Usage: gt cubics.bin out.bin NMAX MODE TOLMIN [J] [M]
 *
 *   cubics.bin  N records of 8 doubles: x0 y0 x1 y1 x2 y2 x3 y3
 *   out.bin     N * NMAX records of NF doubles (n = 1..NMAX), see below
 *   MODE        0: cu2qu ramp  s_k = 1 - 2k/(n-1); n = 1 is the tangent
 *                  intersection, as in fontTools (absent -> +inf)
 *               1: s_k = 0 for every k and n (n = 1: the "midpoint" control)
 *               2: alternating s = +1,-1,+1,... (even n only; n = 1 as mode 0;
 *                  odd n > 1 -> +inf)
 *   TOLMIN      smallest tolerance of interest: n beyond the first one whose
 *               equal-parameter error is <= TOLMIN is not computed (NaN)
 *   J           samples per quadratic for the outer maximum (default 64)
 *   M           ignored (kept for compatibility with gt_scan.c)
 *
 * Fields per (cubic, n):
 *   0 ep     max_t |S(t) - B(t)|                 (equal-parameter error)
 *   1 h_sb   max_t dist(S(t), B)                 (spline -> cubic, directed)
 *   2 h_bs   max_t dist(B(t), S)                 (cubic -> spline, directed)
 *   3 fo_b   max_t |E(t) x B'(t)| / |B'(t)|      (error across the cubic's tangent)
 *   4 fo_s   max_t |E(t) x S'(t)| / |S'(t)|      (error across the spline's tangent)
 *   5 t_sb   parameter where h_sb is attained
 *   6 t_bs   parameter where h_bs is attained
 *   7 t_ep   parameter where ep is attained
 *
 * How distances are computed.  dist(p, C) is exact: the stationary points of
 * |C(tau) - p|^2 are the roots of a quintic (C the cubic) or of one cubic per
 * quadratic (C the spline); all real roots in [0,1] are isolated between the
 * roots of the derivative polynomial (found the same way, recursively) and
 * refined by bisection; the distance is the smallest over those roots and the
 * end points.  The outer maximum over t is taken on J+1 points per quadratic
 * and refined by golden-section search around the three best local maxima,
 * and in nested brackets (down to 1e-7 in t) around the points where the curve
 * the points are taken from has least speed (hairpins and cusps, where the
 * distance can have a peak much narrower than the grid).
 * (gt_scan.c is the first version, which scanned C instead and could settle
 * on a non-global minimum at a hairpin.)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define NF 8

typedef struct { double x, y; } pt;
static inline pt P (double x, double y) { pt p = {x, y}; return p; }
static inline pt add (pt a, pt b) { return P (a.x + b.x, a.y + b.y); }
static inline pt sub (pt a, pt b) { return P (a.x - b.x, a.y - b.y); }
static inline pt mul (pt a, double s) { return P (a.x * s, a.y * s); }
static inline double dot (pt a, pt b) { return a.x * b.x + a.y * b.y; }
static inline double cross (pt a, pt b) { return a.x * b.y - a.y * b.x; }
static inline double norm (pt a) { return hypot (a.x, a.y); }
static inline double dist2 (pt a, pt b) { double dx = a.x - b.x, dy = a.y - b.y; return dx * dx + dy * dy; }

typedef struct { pt c0, c3, b1, b2, b3; } cubic_t;

static inline pt cub_at (const cubic_t *c, double t)
{
  if (t <= 0) return c->c0;
  if (t >= 1) return c->c3;
  return P (c->c0.x + t * (c->b1.x + t * (c->b2.x + t * c->b3.x)),
	    c->c0.y + t * (c->b1.y + t * (c->b2.y + t * c->b3.y)));
}
static inline pt cub_d (const cubic_t *c, double t)
{
  return P (c->b1.x + t * (2 * c->b2.x + 3 * t * c->b3.x),
	    c->b1.y + t * (2 * c->b2.y + 3 * t * c->b3.y));
}

#define MAXN 128
typedef struct { int n; pt on[MAXN + 1]; pt off[MAXN]; } spline_t;

static inline pt spl_at (const spline_t *s, double t)
{
  if (t <= 0) return s->on[0];
  if (t >= 1) return s->on[s->n];
  double x = t * s->n;
  int k = (int) x;
  if (k >= s->n) k = s->n - 1;
  double u = x - k, v = 1 - u;
  pt a = s->on[k], q = s->off[k], b = s->on[k + 1];
  return P (v * v * a.x + 2 * u * v * q.x + u * u * b.x,
	    v * v * a.y + 2 * u * v * q.y + u * u * b.y);
}
static inline pt spl_d (const spline_t *s, double t)
{
  double x = t * s->n;
  int k = (int) x;
  if (k >= s->n) k = s->n - 1;
  if (k < 0) k = 0;
  double u = x - k, v = 1 - u;
  pt a = s->on[k], q = s->off[k], b = s->on[k + 1];
  double f = 2.0 * s->n;
  return P (f * (v * (q.x - a.x) + u * (b.x - q.x)),
	    f * (v * (q.y - a.y) + u * (b.y - q.y)));
}

static inline double horner (const double *c, int deg, double x)
{
  double r = c[deg];
  for (int i = deg - 1; i >= 0; i--) r = r * x + c[i];
  return r;
}

/* All roots in [0,1] of a polynomial of degree <= 5 at which it changes sign
 * (or is exactly zero at an end of a monotone stretch), in increasing order. */
static int roots01 (const double *c, int deg, double *out)
{
  while (deg > 0 && c[deg] == 0) deg--;
  if (deg == 0) return 0;
  if (deg == 1)
  {
    double r = -c[0] / c[1];
    if (r >= 0 && r <= 1) { out[0] = r; return 1; }
    return 0;
  }
  double dc[6], crit[6];
  for (int i = 1; i <= deg; i++) dc[i - 1] = i * c[i];
  int nc = roots01 (dc, deg - 1, crit);
  double pts[8];
  int np = 0;
  pts[np++] = 0;
  for (int i = 0; i < nc; i++) if (crit[i] > 0 && crit[i] < 1) pts[np++] = crit[i];
  pts[np++] = 1;
  int nr = 0;
  for (int i = 0; i + 1 < np; i++)
  {
    double a = pts[i], b = pts[i + 1];
    double fa = horner (c, deg, a), fb = horner (c, deg, b);
    if (fa == 0) { if (!nr || out[nr - 1] != a) out[nr++] = a; continue; }
    if ((fa < 0) == (fb < 0) || fb == 0) continue;
    for (int it = 0; it < 60; it++)
    {
      double m = 0.5 * (a + b);
      double fm = horner (c, deg, m);
      if ((fm < 0) == (fa < 0)) { a = m; fa = fm; } else b = m;
    }
    out[nr++] = 0.5 * (a + b);
  }
  if (horner (c, deg, 1.0) == 0 && (!nr || out[nr - 1] != 1.0)) out[nr++] = 1.0;
  return nr;
}

/* exact distance from p to the cubic */
static double nearest_cubic (const cubic_t *c, pt p, double *tau_out)
{
  pt e0 = sub (c->c0, p);
  double k[6];
  k[0] = dot (e0, c->b1);
  k[1] = 2 * dot (e0, c->b2) + dot (c->b1, c->b1);
  k[2] = 3 * dot (e0, c->b3) + 3 * dot (c->b1, c->b2);
  k[3] = 4 * dot (c->b1, c->b3) + 2 * dot (c->b2, c->b2);
  k[4] = 5 * dot (c->b2, c->b3);
  k[5] = 3 * dot (c->b3, c->b3);
  double r[8];
  int nr = roots01 (k, 5, r);
  double best = dist2 (c->c0, p), tb = 0;
  double d = dist2 (c->c3, p);
  if (d < best) { best = d; tb = 1; }
  for (int i = 0; i < nr; i++)
  {
    d = dist2 (cub_at (c, r[i]), p);
    if (d < best) { best = d; tb = r[i]; }
  }
  if (tau_out) *tau_out = tb;
  return sqrt (best);
}

/* exact distance from p to the spline */
static double nearest_spline (const spline_t *s, pt p, double *tau_out)
{
  double best = dist2 (s->on[0], p), tb = 0;
  for (int k = 0; k < s->n; k++)
  {
    pt a = sub (s->off[k], s->on[k]), z = sub (s->on[k + 1], s->off[k]);
    pt w = sub (z, a);
    pt e0 = sub (s->on[k], p);
    /* cheap rejection: p is farther from the control triangle's box than the best so far */
    {
      double x0 = fmin (s->on[k].x, fmin (s->off[k].x, s->on[k + 1].x));
      double x1 = fmax (s->on[k].x, fmax (s->off[k].x, s->on[k + 1].x));
      double y0 = fmin (s->on[k].y, fmin (s->off[k].y, s->on[k + 1].y));
      double y1 = fmax (s->on[k].y, fmax (s->off[k].y, s->on[k + 1].y));
      double dx = p.x < x0 ? x0 - p.x : p.x > x1 ? p.x - x1 : 0;
      double dy = p.y < y0 ? y0 - p.y : p.y > y1 ? p.y - y1 : 0;
      if (dx * dx + dy * dy > best) continue;
    }
    /* (Q(u) - p) . Q'(u) / 2,  Q(u) - p = e0 + 2 a u + w u^2,  Q'(u) / 2 = a + w u */
    double kk[4];
    kk[0] = dot (e0, a);
    kk[1] = dot (e0, w) + 2 * dot (a, a);
    kk[2] = 3 * dot (a, w);
    kk[3] = dot (w, w);
    double r[6];
    int nr = roots01 (kk, 3, r);
    double d = dist2 (s->on[k + 1], p);
    if (d < best) { best = d; tb = (k + 1.0) / s->n; }
    for (int i = 0; i < nr; i++)
    {
      double u = r[i], v = 1 - u;
      pt q = P (v * v * s->on[k].x + 2 * u * v * s->off[k].x + u * u * s->on[k + 1].x,
		v * v * s->on[k].y + 2 * u * v * s->off[k].y + u * u * s->on[k + 1].y);
      d = dist2 (q, p);
      if (d < best) { best = d; tb = (k + u) / s->n; }
    }
  }
  if (tau_out) *tau_out = tb;
  return sqrt (best);
}

#define GOLD 0.6180339887498949

/* generic "maximise g over [0,1]": samples at i / NS, then golden around the best local maxima */
typedef double (*objfn) (void *ctx, double t);

static double golden_max (objfn g, void *ctx, double lo, double hi, double *targ)
{
  double a = lo, b = hi;
  double x1 = b - GOLD * (b - a), x2 = a + GOLD * (b - a);
  double f1 = g (ctx, x1), f2 = g (ctx, x2);
  for (int it = 0; it < 40; it++)
  {
    if (f1 > f2) { b = x2; x2 = x1; f2 = f1; x1 = b - GOLD * (b - a); f1 = g (ctx, x1); }
    else         { a = x1; x1 = x2; f1 = f2; x2 = a + GOLD * (b - a); f2 = g (ctx, x2); }
  }
  double x = 0.5 * (a + b);
  double f = g (ctx, x);
  *targ = x;
  return f;
}

static double maximise (objfn g, void *ctx, int NS, double *vals, double *targ,
			const double *special, int nspecial)
{
  for (int i = 0; i <= NS; i++) vals[i] = g (ctx, (double) i / NS);
  /* three best local maxima */
  int idx[3] = {-1, -1, -1};
  for (int i = 0; i <= NS; i++)
  {
    double v = vals[i];
    if (!(v == v)) continue;
    if (i > 0 && vals[i - 1] > v) continue;
    if (i < NS && vals[i + 1] > v) continue;
    for (int j = 0; j < 3; j++)
      if (idx[j] < 0 || v > vals[idx[j]])
      {
	for (int l = 2; l > j; l--) idx[l] = idx[l - 1];
	idx[j] = i;
	break;
      }
  }
  double best = -INFINITY, tb = 0;
  for (int j = 0; j < 3; j++)
  {
    int i = idx[j];
    if (i < 0) continue;
    if (vals[i] > best) { best = vals[i]; tb = (double) i / NS; }
    double lo = i > 0 ? (double) (i - 1) / NS : 0.0;
    double hi = i < NS ? (double) (i + 1) / NS : 1.0;
    double tt;
    double f = golden_max (g, ctx, lo, hi, &tt);
    if (f > best) { best = f; tb = tt; }
  }
  /* Special points: where the curve the points are taken from turns around most sharply
   * (speed minima: hairpins, cusps).  The distance can have a peak there far narrower
   * than the grid, so search nested brackets around each. */
  for (int j = 0; j < nspecial; j++)
  {
    double t0 = special[j];
    if (!(t0 >= 0 && t0 <= 1)) continue;
    double v = g (ctx, t0);
    if (v > best) { best = v; tb = t0; }
    for (double h = 1.0 / NS; h > 1e-7; h /= 16)
    {
      double lo = t0 - h > 0 ? t0 - h : 0, hi = t0 + h < 1 ? t0 + h : 1;
      double tt;
      double f = golden_max (g, ctx, lo, hi, &tt);
      if (f > best) { best = f; tb = tt; }
    }
  }
  if (targ) *targ = tb;
  return best;
}

typedef struct { const cubic_t *c; const spline_t *s; } ctx_t;

static double g_ep (void *v, double t)
{ ctx_t *x = (ctx_t *) v; return norm (sub (spl_at (x->s, t), cub_at (x->c, t))); }
static double g_sb (void *v, double t)
{ ctx_t *x = (ctx_t *) v; return nearest_cubic (x->c, spl_at (x->s, t), NULL); }
static double g_bs (void *v, double t)
{ ctx_t *x = (ctx_t *) v; return nearest_spline (x->s, cub_at (x->c, t), NULL); }
static double g_fob (void *v, double t)
{
  ctx_t *x = (ctx_t *) v;
  pt e = sub (spl_at (x->s, t), cub_at (x->c, t));
  pt d = cub_d (x->c, t);
  double nd = norm (d);
  return nd > 0 ? fabs (cross (e, d)) / nd : norm (e);
}
static double g_fos (void *v, double t)
{
  ctx_t *x = (ctx_t *) v;
  pt e = sub (spl_at (x->s, t), cub_at (x->c, t));
  pt d = spl_d (x->s, t);
  double nd = norm (d);
  return nd > 0 ? fabs (cross (e, d)) / nd : norm (e);
}

/* fontTools' calc_intersect + the rule of cubic_approx_quadratic */
static int tangent_intersection (const pt c[4], pt *q)
{
  pt ab = sub (c[1], c[0]), cd = sub (c[3], c[2]);
  pt p = P (-ab.y, ab.x);
  double den = dot (p, cd);
  if (fabs (den) < 1e-15) den = 0;
  if (den == 0)
  {
    int b_eq_c = c[1].x == c[2].x && c[1].y == c[2].y;
    int a_eq_b = c[0].x == c[1].x && c[0].y == c[1].y;
    int c_eq_d = c[2].x == c[3].x && c[2].y == c[3].y;
    if (b_eq_c && (a_eq_b || c_eq_d)) { *q = c[1]; return 1; }
    return 0;
  }
  double num = dot (p, sub (c[0], c[2]));
  if (fabs (num) < 1e-15) num = 0;
  double h = num / den;
  *q = add (c[2], mul (cd, h));
  return isfinite (q->x) && isfinite (q->y);
}

static double weight (int mode, int k, int n)
{
  if (mode == 1) return 0.0;
  if (mode == 2) return (k & 1) ? -1.0 : 1.0;
  return n > 1 ? 1.0 - 2.0 * k / (n - 1) : 0.0;
}

/* build the spline; returns 0 if not available */
static int build_spline (const pt c[4], const cubic_t *cb, int n, int mode, spline_t *s)
{
  s->n = n;
  s->on[0] = c[0];
  s->on[n] = c[3];
  if (n == 1 && mode != 1)
  {
    pt q;
    if (!tangent_intersection (c, &q)) return 0;
    s->off[0] = q;
    return 1;
  }
  if (mode == 2 && (n & 1)) return 0;
  pt d = mul (cb->b3, 1.0 / ((double) n * n * n));
  for (int k = 0; k < n; k++)
  {
    double t0 = (double) k / n, t1 = (double) (k + 1) / n;
    pt p0 = k ? cub_at (cb, t0) : c[0];
    pt p3 = k + 1 < n ? cub_at (cb, t1) : c[3];
    pt d0 = cub_d (cb, t0), d1 = cub_d (cb, t1);
    double third = 1.0 / (3.0 * n);
    pt p1 = add (p0, mul (d0, third));
    pt p2 = sub (p3, mul (d1, third));
    /* q_M = (3 (p1 + p2) - p0 - p3) / 4 */
    pt qm = P ((3 * (p1.x + p2.x) - p0.x - p3.x) * 0.25, (3 * (p1.y + p2.y) - p0.y - p3.y) * 0.25);
    double w = weight (mode, k, n);
    s->off[k] = add (qm, mul (d, w * 0.25));
  }
  for (int k = 1; k < n; k++)
    s->on[k] = mul (add (s->off[k - 1], s->off[k]), 0.5);
  return 1;
}

int main (int argc, char **argv)
{
  if (argc < 6) { fprintf (stderr, "usage: gt cubics.bin out.bin NMAX MODE TOLMIN [J] [M]\n"); return 1; }
  int NMAX = atoi (argv[3]);
  int mode = atoi (argv[4]);
  double tolmin = atof (argv[5]);
  int J = argc > 6 ? atoi (argv[6]) : 64;
  int M = argc > 7 ? atoi (argv[7]) : 256;
  if (NMAX > MAXN) NMAX = MAXN;
  (void) M; /* unused since the nearest-point search became exact */

  FILE *f = fopen (argv[1], "rb");
  if (!f) { perror (argv[1]); return 1; }
  fseek (f, 0, SEEK_END);
  long sz = ftell (f);
  fseek (f, 0, SEEK_SET);
  long N = sz / (8 * sizeof (double));
  double *in = malloc (sz);
  if (fread (in, 8 * sizeof (double), N, f) != (size_t) N) { fprintf (stderr, "short read\n"); return 1; }
  fclose (f);

  double *out = malloc ((size_t) N * NMAX * NF * sizeof (double));
  for (long i = 0; i < N * NMAX * NF; i++) out[i] = NAN;

#pragma omp parallel
  {
    double *vals = malloc (((size_t) J * NMAX + 1) * sizeof (double));
    spline_t *s = malloc (sizeof (spline_t));

#pragma omp for schedule(dynamic, 64)
    for (long i = 0; i < N; i++)
    {
      pt c[4];
      for (int j = 0; j < 4; j++) c[j] = P (in[i * 8 + 2 * j], in[i * 8 + 2 * j + 1]);
      cubic_t cb;
      cb.c0 = c[0]; cb.c3 = c[3];
      cb.b1 = mul (sub (c[1], c[0]), 3);
      cb.b2 = P (3 * (c[0].x - 2 * c[1].x + c[2].x), 3 * (c[0].y - 2 * c[1].y + c[2].y));
      cb.b3 = P (c[3].x - 3 * c[2].x + 3 * c[1].x - c[0].x, c[3].y - 3 * c[2].y + 3 * c[1].y - c[0].y);

      /* speed minima of the cubic: roots of B'.B'' */
      double sp_b[4];
      int ns_b;
      {
	double kk[4];
	kk[0] = 2 * dot (cb.b1, cb.b2);
	kk[1] = 6 * dot (cb.b1, cb.b3) + 4 * dot (cb.b2, cb.b2);
	kk[2] = 18 * dot (cb.b2, cb.b3);
	kk[3] = 18 * dot (cb.b3, cb.b3);
	ns_b = roots01 (kk, 3, sp_b);
      }

      for (int n = 1; n <= NMAX; n++)
      {
	double *o = out + ((size_t) i * NMAX + (n - 1)) * NF;
	if (!build_spline (c, &cb, n, mode, s))
	{
	  for (int j = 0; j < 5; j++) o[j] = INFINITY;
	  continue;
	}
	int NS = J * n;
	ctx_t x_sb = { &cb, s };
	ctx_t x_bs = { &cb, s };

	/* speed minima of the spline: the vertex of each quadratic */
	double sp_s[MAXN];
	int ns_s = 0;
	for (int k = 0; k < n; k++)
	{
	  pt a = sub (s->off[k], s->on[k]), w = sub (sub (s->on[k + 1], s->off[k]), a);
	  double ww = dot (w, w);
	  if (ww > 0)
	  {
	    double u = -dot (a, w) / ww;
	    if (u > 0 && u < 1) sp_s[ns_s++] = (k + u) / n;
	  }
	}

	o[0] = maximise (g_ep, &x_sb, NS, vals, &o[7], NULL, 0);
	o[1] = maximise (g_sb, &x_sb, NS, vals, &o[5], sp_s, ns_s);
	o[2] = maximise (g_bs, &x_bs, NS, vals, &o[6], sp_b, ns_b);
	o[3] = maximise (g_fob, &x_sb, NS, vals, NULL, NULL, 0);
	o[4] = maximise (g_fos, &x_sb, NS, vals, NULL, NULL, 0);

	if (o[0] <= tolmin) break;
      }
    }
    free (vals); free (s);
  }

  f = fopen (argv[2], "wb");
  if (!f) { perror (argv[2]); return 1; }
  fwrite (out, sizeof (double), (size_t) N * NMAX * NF, f);
  fclose (f);
  fprintf (stderr, "%ld cubics, NMAX %d, mode %d -> %s\n", N, NMAX, mode, argv[2]);
  return 0;
}
