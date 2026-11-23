/* Hausdorff distance between a cubic Bezier and an implied-on-curve quadratic
 * spline: a measuring instrument (dense sampling, Newton for the nearest
 * point, golden-section refinement of each local maximum).
 *
 * Build:  gcc -O2 -fopenmp -shared -fPIC -o libchd.so chd.c -lm
 *
 * Cubic c: 8 doubles (x0 y0 .. x3 y3).  Spline: n off-curve points Q (2n
 * doubles); quadratic k is on_k, Q_k, on_{k+1} with on_0 = c0, on_n = c3,
 * on_k = (Q_{k-1} + Q_k) / 2.  Spline parameter xi in [0, n].
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

static inline void cubic_eval (const double *c, double t, double *P, double *D, double *DD)
{
  for (int a = 0; a < 2; a++)
  {
    double c0 = c[a], c1 = c[2 + a], c2 = c[4 + a], c3 = c[6 + a];
    double b1 = 3 * (c1 - c0), b2 = 3 * (c0 - 2 * c1 + c2), b3 = c3 - 3 * c2 + 3 * c1 - c0;
    P[a] = c0 + t * (b1 + t * (b2 + t * b3));
    if (D) D[a] = b1 + t * (2 * b2 + 3 * t * b3);
    if (DD) DD[a] = 2 * b2 + 6 * t * b3;
  }
}

static inline void spline_eval (const double *c, const double *Q, int n, double xi,
				double *P, double *D, double *DD)
{
  int k = (int) floor (xi);
  if (k < 0) k = 0;
  if (k > n - 1) k = n - 1;
  double u = xi - k;
  for (int a = 0; a < 2; a++)
  {
    double q = Q[2 * k + a];
    double s = k == 0 ? c[a] : 0.5 * (Q[2 * (k - 1) + a] + q);
    double e = k == n - 1 ? c[6 + a] : 0.5 * (q + Q[2 * (k + 1) + a]);
    P[a] = (1 - u) * (1 - u) * s + 2 * u * (1 - u) * q + u * u * e;
    if (D) D[a] = 2 * (1 - u) * (q - s) + 2 * u * (e - q);
    if (DD) DD[a] = 2 * (s - 2 * q + e);
  }
}

/* Nearest point of the cubic to p, searching [lo, hi] from t. */
static double newton_cubic (const double *c, const double *p, double t, double lo, double hi, double *t_out)
{
  double P[2], D[2], DD[2];
  for (int it = 0; it < 8; it++)
  {
    cubic_eval (c, t, P, D, DD);
    double vx = P[0] - p[0], vy = P[1] - p[1];
    double f = vx * D[0] + vy * D[1];
    double fp = D[0] * D[0] + D[1] * D[1] + vx * DD[0] + vy * DD[1];
    double step = fp > 1e-30 ? f / fp : (f > 0 ? 1 : -1) * (hi - lo) * 0.25;
    double tn = t - step;
    if (tn < lo) tn = lo;
    if (tn > hi) tn = hi;
    if (fabs (tn - t) < 1e-14) { t = tn; break; }
    t = tn;
  }
  cubic_eval (c, t, P, NULL, NULL);
  *t_out = t;
  return hypot (P[0] - p[0], P[1] - p[1]);
}

static double nearest_cubic (const double *c, const double *Bs, int K, const double *p, double *t_out)
{
  int j = 0;
  double best = 1e300;
  for (int i = 0; i < K; i++)
  {
    double dx = Bs[2 * i] - p[0], dy = Bs[2 * i + 1] - p[1];
    double d = dx * dx + dy * dy;
    if (d < best) { best = d; j = i; }
  }
  double h = 1.0 / (K - 1), t0 = j * h;
  double lo = t0 - h < 0 ? 0 : t0 - h, hi = t0 + h > 1 ? 1 : t0 + h;
  double t;
  double d = newton_cubic (c, p, t0, lo, hi, &t);
  if (d > sqrt (best)) { d = sqrt (best); t = t0; }
  *t_out = t;
  return d;
}

static double newton_spline (const double *c, const double *Q, int n, const double *p,
			     double x, double lo, double hi, double *x_out)
{
  double P[2], D[2], DD[2];
  for (int it = 0; it < 8; it++)
  {
    spline_eval (c, Q, n, x, P, D, DD);
    double vx = P[0] - p[0], vy = P[1] - p[1];
    double f = vx * D[0] + vy * D[1];
    double fp = D[0] * D[0] + D[1] * D[1] + vx * DD[0] + vy * DD[1];
    double step = fp > 1e-30 ? f / fp : (f > 0 ? 1 : -1) * (hi - lo) * 0.25;
    double xn = x - step;
    if (xn < lo) xn = lo;
    if (xn > hi) xn = hi;
    if (fabs (xn - x) < 1e-13) { x = xn; break; }
    x = xn;
  }
  spline_eval (c, Q, n, x, P, NULL, NULL);
  *x_out = x;
  return hypot (P[0] - p[0], P[1] - p[1]);
}

static double nearest_spline (const double *c, const double *Q, int n, const double *Ss, int K,
			      const double *p, double *x_out)
{
  int j = 0;
  double best = 1e300;
  for (int i = 0; i < K; i++)
  {
    double dx = Ss[2 * i] - p[0], dy = Ss[2 * i + 1] - p[1];
    double d = dx * dx + dy * dy;
    if (d < best) { best = d; j = i; }
  }
  double h = (double) n / (K - 1), x0 = j * h;
  double lo = x0 - h < 0 ? 0 : x0 - h, hi = x0 + h > n ? n : x0 + h;
  double x;
  double d = newton_spline (c, Q, n, p, x0, lo, hi, &x);
  if (d > sqrt (best)) { d = sqrt (best); x = x0; }
  *x_out = x;
  return d;
}

/* Golden-section maximum of f on [lo, hi]; f is the distance from a point of
 * one curve (at the given parameter) to the other curve. */
typedef struct {
  const double *c, *Q; int n;
  const double *Bs; int K;      /* coarse table of the cubic */
  const double *Ss; int Ks;     /* coarse table of the spline */
} ctx_t;

static double dist_s2c (const ctx_t *x, double xi)
{
  double P[2], t;
  spline_eval (x->c, x->Q, x->n, xi, P, NULL, NULL);
  return nearest_cubic (x->c, x->Bs, x->K, P, &t);
}
static double dist_c2s (const ctx_t *x, double t)
{
  double P[2], xi;
  cubic_eval (x->c, t, P, NULL, NULL);
  return nearest_spline (x->c, x->Q, x->n, x->Ss, x->Ks, P, &xi);
}

static double golden_max (const ctx_t *x, double (*f) (const ctx_t *, double),
			  double lo, double hi, double *where)
{
  const double gr = 0.6180339887498949;
  double x1 = hi - gr * (hi - lo), x2 = lo + gr * (hi - lo);
  double f1 = f (x, x1), f2 = f (x, x2);
  for (int it = 0; it < 16; it++)
  {
    if (f1 < f2) { lo = x1; x1 = x2; f1 = f2; x2 = lo + gr * (hi - lo); f2 = f (x, x2); }
    else { hi = x2; x2 = x1; f2 = f1; x1 = hi - gr * (hi - lo); f1 = f (x, x1); }
  }
  if (where) *where = f1 > f2 ? x1 : x2;
  return f1 > f2 ? f1 : f2;
}

/* max over [0, len] of f, from M = steps + 1 equally spaced samples: every
 * local maximum among the samples is refined, and so are the first and the
 * last interval (both curves share their end points, so the samples there are
 * zero and can hide a bump). */
static double scan_max (const ctx_t *x, double (*f) (const ctx_t *, double),
			double len, int steps, double *d, double *where)
{
  double h = len / steps, best = 0, w = 0;
  for (int i = 0; i <= steps; i++) { d[i] = f (x, i * h); if (d[i] > best) { best = d[i]; w = i * h; } }
  for (int i = 0; i <= steps; i++)
  {
    double lo, hi;
    if (i == 0) { lo = 0; hi = h; }
    else if (i == steps) { lo = len - h; hi = len; }
    else if (d[i] >= d[i - 1] && d[i] >= d[i + 1] && d[i] > 0.5 * best) { lo = (i - 1) * h; hi = (i + 1) * h; }
    else continue;
    double ww, v = golden_max (x, f, lo, hi, &ww);
    if (v > best) { best = v; w = ww; }
  }
  if (where) *where = w;
  return best;
}

/* Directed distances and Hausdorff distance.  m samples per quadratic. */
double hd_one (const double *c, const double *Q, int n, int m, int both,
	       double *s2c, double *c2s, double *where)
{
  int steps = n * m;
  /* coarse tables for the nearest-point search: at least 64 intervals, and
   * twice as fine as the evaluation samples */
  int K = 2 * steps; if (K < 64) K = 64; K += 1;
  int sm = 2 * m; if (n * sm < 64) sm = (64 + n - 1) / n;
  int Ks = n * sm + 1;
  double *Bs = (double *) malloc (sizeof (double) * 2 * (size_t) K);
  double *Ss = (double *) malloc (sizeof (double) * 2 * (size_t) Ks);
  double *d = (double *) malloc (sizeof (double) * (size_t) (steps + 1));
  for (int i = 0; i < K; i++) cubic_eval (c, (double) i / (K - 1), Bs + 2 * i, NULL, NULL);
  for (int i = 0; i < Ks; i++) spline_eval (c, Q, n, (double) i / sm, Ss + 2 * i, NULL, NULL);
  ctx_t x = { c, Q, n, Bs, K, Ss, Ks };

  double awhere = 0;
  double a = scan_max (&x, dist_s2c, (double) n, steps, d, &awhere);
  double b = both ? scan_max (&x, dist_c2s, 1.0, steps, d, NULL) : 0;
  free (Bs); free (Ss); free (d);
  if (s2c) *s2c = a;
  if (c2s) *c2s = b;
  if (where) *where = awhere;
  return a > b ? a : b;
}

/* N cubics (8 doubles each) and N splines of n off-curve points each. */
void hd_batch (const double *C, const double *Q, long N, int n, int m, int both, double *out)
{
#pragma omp parallel for schedule(dynamic, 64)
  for (long i = 0; i < N; i++)
    out[i] = hd_one (C + 8 * i, Q + 2 * (size_t) n * i, n, m, both, NULL, NULL, NULL);
}

/* Largest distance at equal parameters: |S(k + u) - B((k + u) / n)|. */
void ep_batch (const double *C, const double *Q, long N, int n, int m, double *out)
{
#pragma omp parallel for schedule(dynamic, 256)
  for (long i = 0; i < N; i++)
  {
    const double *c = C + 8 * i, *q = Q + 2 * (size_t) n * i;
    double best = 0;
    for (int j = 0; j <= n * m; j++)
    {
      double x = (double) j / m, P[2], B[2];
      spline_eval (c, q, n, x, P, NULL, NULL);
      cubic_eval (c, x / n, B, NULL, NULL);
      double d = hypot (P[0] - B[0], P[1] - B[1]);
      if (d > best) best = d;
    }
    out[i] = best;
  }
}

/* Residuals for an optimiser: at xi_i = i / m, i = 0..n m, the signed distance
 * r_i from the spline point to the cubic and its gradient (gx, gy) with respect
 * to the spline point.  Interior foot: r = (P - B(t)) . N(t), gradient N(t).
 * Foot at an end of the cubic: r = |P - end|, gradient the unit vector. */
void resid (const double *c, const double *Q, int n, int m, double *r, double *gx, double *gy, double *tf)
{
  int M = n * m + 1;
  int K = 2 * n * m + 1;
  if (K < 65) K = 65;
  double *Bs = (double *) malloc (sizeof (double) * 2 * (size_t) K);
  for (int i = 0; i < K; i++) cubic_eval (c, (double) i / (K - 1), Bs + 2 * i, NULL, NULL);
  for (int i = 0; i < M; i++)
  {
    double P[2], B[2], D[2], t;
    spline_eval (c, Q, n, (double) i / m, P, NULL, NULL);
    double d = nearest_cubic (c, Bs, K, P, &t);
    cubic_eval (c, t, B, D, NULL);
    double len = hypot (D[0], D[1]);
    double vx = P[0] - B[0], vy = P[1] - B[1];
    if (t > 0 && t < 1 && len > 1e-12)
    {
      double nx = -D[1] / len, ny = D[0] / len;
      r[i] = vx * nx + vy * ny; gx[i] = nx; gy[i] = ny;
    }
    else if (d > 1e-12)
    {
      r[i] = d; gx[i] = vx / d; gy[i] = vy / d;
    }
    else
    {
      double nx = len > 1e-12 ? -D[1] / len : 0, ny = len > 1e-12 ? D[0] / len : 0;
      r[i] = 0; gx[i] = nx; gy[i] = ny;
    }
    if (tf) tf[i] = t;
  }
  free (Bs);
}
