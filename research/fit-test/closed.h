/* Closed-form versions of the fit test for a single quadratic: is the
 * cubic Bezier (0, e1, e2, 0) within tol of the origin?
 *
 * With m = (e1 + e2)/2, r = (e2 - e1)/2 and t = 2u - 1 in [-1, 1] the curve
 * is (3/4) (1 - t^2) (m + t r). With A = |m|^2, B = m.r, C = |r|^2 its
 * squared distance is (9/16) G(t), G(t) = (1 - t^2)^2 (A + 2 B t + C t^2), and
 * G'(t) = -2 (1 - t^2) k(t), k(t) = 3C t^3 + 5B t^2 + (2A - C) t - B.
 *
 * Take B >= 0 (otherwise look from the other end: t -> -t). Then
 * G(t) >= G(-t) for t >= 0, so the largest distance is in [0, 1]; and
 * k(0) = -B <= 0, k(1) = 2 |e2|^2 >= 0, k'' = 18C t + 10B >= 0: k is convex
 * on [0, 1] and has one root there (the larger one if B = 0), where G is
 * largest. Newton's iteration from t = 1 comes down to it and never passes it.
 *
 * Only + - * / and comparisons, in a fixed order: the same everywhere.
 */
#ifndef CLOSED_H
#define CLOSED_H

#include "fitsim.h"

static int closed_iters; /* Newton steps of the last call */

/* Stop when a Newton step can no longer raise G by a relative 2^-55. */
#define CLOSED_MAX_STEPS 64

static inline double closed_root(double A, double B, double C)
{
    /* B >= 0 */
    const double c3 = 3.0 * C, c2 = 5.0 * B, c1 = 2.0 * A - C;
    const double d2 = 9.0 * C, d1 = 10.0 * B;
    const double small = (A + C) * 0x1p-56;
    double t = 1.0;
    int i;
    for (i = 0; i < CLOSED_MAX_STEPS; i++) {
        double k = ((c3 * t + c2) * t + c1) * t - B;
        if (!(k > 0.0)) break;
        double dk = (d2 * t + d1) * t + c1;
        double step = k / dk;
        t -= step;
        if (!(k * step > small)) break;
    }
    closed_iters = i;
    if (!(t > 0.0)) t = 0.0;
    return t;
}

/* The whole test, with the two quick answers of the recursive test first
 * (they are its first level, operation for operation). */
static int closed_fit(V e1, V e2, double tol)
{
    closed_iters = -1;
    if (vabs(e2) <= tol && vabs(e1) <= tol) return 1;
    V zero = { 0.0, 0.0 };
    V mid = vmul(vadd(vadd(zero, vmul(vadd(e1, e2), 3.0)), zero), 0.125);
    if (!(vabs(mid) <= tol)) return 0;
    V m = vmul(vadd(e1, e2), 0.5), r = vmul(vsub(e2, e1), 0.5);
    double A = m.x * m.x + m.y * m.y;
    double B = m.x * r.x + m.y * r.y;
    double C = r.x * r.x + r.y * r.y;
    if (B < 0.0) { B = -B; r.x = -r.x; r.y = -r.y; }
    double t = closed_root(A, B, C);
    V v = vadd(m, vmul(r, t));
    return 0.75 * (1.0 - t * t) * vabs(v) <= tol;
}

/* Without the quick answers: the closed form alone. */
static int closed_only(V e1, V e2, double tol)
{
    V m = vmul(vadd(e1, e2), 0.5), r = vmul(vsub(e2, e1), 0.5);
    double A = m.x * m.x + m.y * m.y;
    double B = m.x * r.x + m.y * r.y;
    double C = r.x * r.x + r.y * r.y;
    if (B < 0.0) { B = -B; r.x = -r.x; r.y = -r.y; }
    double t = closed_root(A, B, C);
    V v = vadd(m, vmul(r, t));
    return 0.75 * (1.0 - t * t) * vabs(v) <= tol;
}

/* The largest distance as the closed form computes it, in double. */
static double closed_max(V e1, V e2)
{
    V m = vmul(vadd(e1, e2), 0.5), r = vmul(vsub(e2, e1), 0.5);
    double A = m.x * m.x + m.y * m.y;
    double B = m.x * r.x + m.y * r.y;
    double C = r.x * r.x + r.y * r.y;
    if (B < 0.0) { B = -B; r.x = -r.x; r.y = -r.y; }
    double t = closed_root(A, B, C);
    V v = vadd(m, vmul(r, t));
    return 0.75 * (1.0 - t * t) * vabs(v);
}

#endif
