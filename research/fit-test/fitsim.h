/* A C port of the arithmetic of fontTools' cu2qu fit test for a single
 * quadratic, operation for operation as pure Python does it on complex
 * numbers (and as the Cython module does it): every operation is one IEEE
 * double operation per coordinate, in Python's order; abs() of a complex is
 * libm hypot().
 *
 * Compile with -ffp-contract=off (no FMA), no -ffast-math.
 */
#ifndef FITSIM_H
#define FITSIM_H

#include <math.h>
#include <stdint.h>
#include <string.h>

typedef struct { double x, y; } V;

static inline V vadd(V a, V b) { V r = { a.x + b.x, a.y + b.y }; return r; }
static inline V vsub(V a, V b) { V r = { a.x - b.x, a.y - b.y }; return r; }
static inline V vmul(V a, double s) { V r = { a.x * s, a.y * s }; return r; }
/* -DSQRT_HYPOT: as kurbo's Vec2::hypot, the square root of x*x + y*y. */
#ifdef SQRT_HYPOT
static inline double vabs(V a) { return sqrt(a.x * a.x + a.y * a.y); }
#else
static inline double vabs(V a) { return hypot(a.x, a.y); }
#endif

/* Statistics of the last fit() call. */
static long fit_nodes;
static int fit_depth;
/* The piece at which the depth limit was hit, if it was. */
static V fit_trap[4];

#define FIT_TRUE 1
#define FIT_FALSE 0
#define FIT_DEEP (-1)

/* cubic_farthest_fit_inside, to at most `limit` halvings. */
static int fit_rec(V p0, V p1, V p2, V p3, double tol, int depth, int limit)
{
    fit_nodes++;
    if (depth > fit_depth) fit_depth = depth;
    if (vabs(p2) <= tol && vabs(p1) <= tol) return FIT_TRUE;
    /* mid = (p0 + 3 * (p1 + p2) + p3) * 0.125 */
    V mid = vmul(vadd(vadd(p0, vmul(vadd(p1, p2), 3.0)), p3), 0.125);
    if (!(vabs(mid) <= tol)) return FIT_FALSE;
    if (depth >= limit) {
        fit_trap[0] = p0; fit_trap[1] = p1; fit_trap[2] = p2; fit_trap[3] = p3;
        return FIT_DEEP;
    }
    /* deriv3 = (p3 + p2 - p1 - p0) * 0.125 */
    V deriv3 = vmul(vsub(vsub(vadd(p3, p2), p1), p0), 0.125);
    int r = fit_rec(p0, vmul(vadd(p0, p1), 0.5), vsub(mid, deriv3), mid, tol, depth + 1, limit);
    if (r != FIT_TRUE) return r;
    return fit_rec(mid, vadd(mid, deriv3), vmul(vadd(p2, p3), 0.5), p3, tol, depth + 1, limit);
}

static int fit(V e1, V e2, double tol, int limit)
{
    V zero = { 0.0, 0.0 };
    fit_nodes = 0;
    fit_depth = 0;
    return fit_rec(zero, e1, e2, zero, tol, 0, limit);
}

/* The same with a depth limit after which the answer is true (kurbo's rule),
 * counting nodes; gives up (FIT_DEEP) after max_nodes. */
static int fit_cap_rec(V p0, V p1, V p2, V p3, double tol, int depth, int limit, long max_nodes)
{
    fit_nodes++;
    if (depth > fit_depth) fit_depth = depth;
    if (fit_nodes > max_nodes) return FIT_DEEP;
    if (vabs(p2) <= tol && vabs(p1) <= tol) return FIT_TRUE;
    V mid = vmul(vadd(vadd(p0, vmul(vadd(p1, p2), 3.0)), p3), 0.125);
    if (!(vabs(mid) <= tol)) return FIT_FALSE;
    if (depth == limit) return FIT_TRUE;
    V deriv3 = vmul(vsub(vsub(vadd(p3, p2), p1), p0), 0.125);
    int r = fit_cap_rec(p0, vmul(vadd(p0, p1), 0.5), vsub(mid, deriv3), mid, tol, depth + 1, limit, max_nodes);
    if (r != FIT_TRUE) return r;
    return fit_cap_rec(mid, vadd(mid, deriv3), vmul(vadd(p2, p3), 0.5), p3, tol, depth + 1, limit, max_nodes);
}

static int fit_cap(V e1, V e2, double tol, int limit, long max_nodes)
{
    V zero = { 0.0, 0.0 };
    fit_nodes = 0;
    fit_depth = 0;
    return fit_cap_rec(zero, e1, e2, zero, tol, 0, limit, max_nodes);
}

/* dot() of cu2qu.py: the real part of v1 * conj(v2), zeroed if tiny. */
static inline double cu2qu_dot(V v1, V v2)
{
    double result = v1.x * v2.x - v1.y * (-v2.y);
    if (fabs(result) < 1e-15) result = 0.0;
    return result;
}

/* cubic_approx_quadratic up to the fit test: returns 1 and sets e1, e2 and
 * q if the single quadratic is tried, 0 if it is not (third difference,
 * no intersection). */
static int single_errors(const V c[4], double tol, V *e1, V *e2, V *q)
{
    /* abs(p3 - 3 * p2 + 3 * p1 - p0) */
    V d3 = vsub(vadd(vsub(c[3], vmul(c[2], 3.0)), vmul(c[1], 3.0)), c[0]);
    if (!(vabs(d3) <= tol * 21.0)) return 0;
    /* calc_intersect(a, b, c, d) */
    V ab = vsub(c[1], c[0]);
    V cd = vsub(c[3], c[2]);
    /* p = ab * 1j: (ab.x*0 - ab.y*1, ab.x*1 + ab.y*0) */
    V p = { ab.x * 0.0 - ab.y * 1.0, ab.x * 1.0 + ab.y * 0.0 };
    double den = cu2qu_dot(p, cd);
    V q1;
    if (den == 0.0) {
        if (c[1].x == c[2].x && c[1].y == c[2].y &&
            ((c[0].x == c[1].x && c[0].y == c[1].y) || (c[2].x == c[3].x && c[2].y == c[3].y)))
            q1 = c[1];
        else
            return 0;
    } else {
        double h = cu2qu_dot(p, vsub(c[0], c[2])) / den;
        q1 = vadd(c[2], vmul(cd, h));
    }
    if (isnan(q1.y)) return 0;
    const double two_thirds = 2.0 / 3.0;
    V c1 = vadd(c[0], vmul(vsub(q1, c[0]), two_thirds));
    V c2 = vadd(c[3], vmul(vsub(q1, c[3]), two_thirds));
    *e1 = vsub(c1, c[1]);
    *e2 = vsub(c2, c[2]);
    *q = q1;
    return 1;
}

/* The largest |3u(1-u)^2 e1 + 3u^2(1-u) e2| over [0, 1], in long double:
 * with m = (e1 + e2)/2, r = (e2 - e1)/2 and t = 2u - 1 it is
 * (3/4) (1 - t^2) |m + t r|, and with A = |m|^2, B = m.r >= 0 (swap the ends
 * otherwise), C = |r|^2 the maximum is at the one root in [0, 1] of the
 * cubic 3C t^3 + 5B t^2 + (2A - C) t - B, which is convex there: Newton
 * from 1 goes down to it. */
static long double max_error_l(V e1, V e2, long double *t_out)
{
    long double mx = ((long double)e1.x + e2.x) / 2, my = ((long double)e1.y + e2.y) / 2;
    long double rx = ((long double)e2.x - e1.x) / 2, ry = ((long double)e2.y - e1.y) / 2;
    long double A = mx * mx + my * my, B = mx * rx + my * ry, C = rx * rx + ry * ry;
    if (B < 0) { B = -B; rx = -rx; ry = -ry; }
    long double t = 1.0L;
    if (A + C == 0) { if (t_out) *t_out = 0; return 0; }
    for (int i = 0; i < 200; i++) {
        long double k = ((3 * C * t + 5 * B) * t + (2 * A - C)) * t - B;
        long double dk = (9 * C * t + 10 * B) * t + (2 * A - C);
        if (!(dk > 0)) break;
        long double tn = t - k / dk;
        if (!(tn < t)) break;
        t = tn;
        if (t < 0) { t = 0; break; }
    }
    long double vx = mx + t * rx, vy = my + t * ry;
    long double best = 0.75L * (1 - t * t) * hypotl(vx, vy);
    if (t_out) *t_out = t;
    return best;
}

static inline uint64_t dbits(double x) { uint64_t u; memcpy(&u, &x, 8); return u; }
static inline double dfrom(uint64_t u) { double x; memcpy(&x, &u, 8); return x; }

/* xorshift64* random numbers */
static inline uint64_t rnd64(uint64_t *s)
{
    uint64_t x = *s;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    *s = x;
    return x * 0x2545F4914F6CDD1DULL;
}
static inline double rnd_unit(uint64_t *s) { return (rnd64(s) >> 11) * (1.0 / 9007199254740992.0); }
static inline long rnd_int(uint64_t *s, long lo, long hi) { return lo + (long)(rnd64(s) % (uint64_t)(hi - lo + 1)); }

#endif
