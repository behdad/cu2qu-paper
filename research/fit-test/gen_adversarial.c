/* Inputs made to sit at the limit of the single quadratic, of the kinds on
 * which the recursive fit test is at its worst: the tolerance is the double
 * nearest the largest distance between the cubic and its single quadratic,
 * and the three doubles on either side.
 *
 * Usage: gen_adversarial <family> <seed> <count>   (writes kdiff_py.py's
 *        format to stdout: "S <tol> <8 coords>", floats as hex bits)
 *
 *   three     three equal points: (P, P, P, D) or (D, P, P, P); the error
 *             curve is on a line through the origin
 *   threeint  the same with integer coordinates
 *   handle    one handle exactly where the single quadratic has it (2/3 of
 *             the way to the tangents' intersection), integer coordinates:
 *             again an error curve on a line
 *   axis      e1 = (p, 1) s, e2 = (0, -2) s: both handles off, the largest
 *             distance on a coordinate axis, at t = 1/3
 *   flat      symmetric cubics whose error is flat to fourth order in the
 *             middle: handles at a slope of 1/sqrt(2) to the chord
 *   ess       e2 = -e1: an S, equally far twice
 *   near      nearly quadratic integer cubics (both handles a little off)
 */
#include <stdio.h>
#include <stdlib.h>
#include "fitsim.h"

static void emit(const V c[4], double tol)
{
    printf("S %016llx", (unsigned long long)dbits(tol));
    for (int i = 0; i < 4; i++)
        printf(" %016llx %016llx", (unsigned long long)dbits(c[i].x), (unsigned long long)dbits(c[i].y));
    printf("\n");
}

int main(int argc, char **argv)
{
    const char *family = argv[1];
    uint64_t seed = strtoull(argv[2], 0, 10) * 0x9E3779B97F4A7C15ULL + 4711;
    long count = atol(argv[3]);
    for (int i = 0; i < 20; i++) rnd64(&seed);
    long made = 0, attempts = 0;
    while (made < count && attempts < count * 50) {
        attempts++;
        V c[4];
        if (!strcmp(family, "three") || !strcmp(family, "threeint")) {
            V P, D;
            if (!strcmp(family, "threeint")) {
                P.x = (double)rnd_int(&seed, -1000, 1000); P.y = (double)rnd_int(&seed, -1000, 1000);
                D.x = P.x + (double)rnd_int(&seed, -300, 300);
                D.y = (rnd64(&seed) & 1) ? P.y : P.y + (double)rnd_int(&seed, -300, 300);
                if (rnd64(&seed) % 4 == 0) { P.x = P.y = 0; }
            } else {
                double s = ldexp(1.0, (int)rnd_int(&seed, 0, 10));
                int at_origin = rnd64(&seed) & 1;
                P.x = at_origin ? 0 : (rnd_unit(&seed) * 2 - 1) * s; P.y = at_origin ? 0 : (rnd_unit(&seed) * 2 - 1) * s;
                D.x = P.x + (rnd_unit(&seed) * 2 - 1) * s;
                D.y = (rnd64(&seed) & 1) ? P.y : P.y + (rnd_unit(&seed) * 2 - 1) * s;
            }
            if (rnd64(&seed) & 1) { c[0] = D; c[1] = c[2] = c[3] = P; }
            else { c[0] = c[1] = c[2] = P; c[3] = D; }
        } else if (!strcmp(family, "handle")) {
            /* q integer; c0 = q - 3a, c3 = q - 3b (a, b integer vectors);
             * one handle at 2/3 (c = end + 2 a or 2 b), the other at j/3 */
            V q = { (double)rnd_int(&seed, -500, 500), (double)rnd_int(&seed, -500, 500) };
            V a = { (double)rnd_int(&seed, -60, 60), (double)rnd_int(&seed, -60, 60) };
            V b = { (double)rnd_int(&seed, -60, 60), (double)rnd_int(&seed, -60, 60) };
            if (rnd64(&seed) % 3 == 0) a.y = 0;
            if (rnd64(&seed) % 3 == 0) b.x = 0;
            if (a.x * b.y - a.y * b.x == 0) continue;
            int j = (int)rnd_int(&seed, 0, 3); /* the other handle at j/3: 0, 1/3, (2/3 is exact), 1 */
            if (j == 2) j = 3;
            c[0] = vsub(q, vmul(a, 3.0));
            c[3] = vsub(q, vmul(b, 3.0));
            c[1] = vadd(c[0], vmul(a, 2.0));
            c[2] = vadd(c[3], vmul(b, 2.0));
            if (rnd64(&seed) & 1) c[1] = vadd(c[0], vmul(a, (double)j));
            else c[2] = vadd(c[3], vmul(b, (double)j));
        } else if (!strcmp(family, "axis")) {
            /* tangents meet at q; q - c0 along (p, 1), q - c3 along (0, 1):
             * c0 = q - 3 (p, 1) u, c1 = c0 + 2 (p,1) u - (p,1) u = q - 2 (p,1) u
             * gives e1 = (p, 1) u; c3 = q + (0, 9) u, c2 = q + (0, 5) u gives
             * e2 = (0, -2) u. Integer p and u, sometimes swapped axes. */
            double p = (double)rnd_int(&seed, 2, 400), u = (double)rnd_int(&seed, 1, 4);
            V q = { (double)rnd_int(&seed, -300, 300), (double)rnd_int(&seed, -300, 300) };
            if (rnd64(&seed) % 3 == 0) { q.x = q.y = 0; }
            V d = { p * u, u };
            c[0] = vsub(q, vmul(d, 3.0));
            c[1] = vsub(q, vmul(d, 2.0));
            c[2].x = q.x; c[2].y = q.y + 5 * u;
            c[3].x = q.x; c[3].y = q.y + 9 * u;
            if (rnd64(&seed) & 1)
                for (int i = 0; i < 4; i++) { double t = c[i].x; c[i].x = c[i].y; c[i].y = t; }
        } else if (!strcmp(family, "flat")) {
            /* c0 = (-X, 0), c3 = (X, 0), handles (X - x1) along and
             * (X - x1)/sqrt(2) up */
            double X = ldexp(1.0 + rnd_unit(&seed), (int)rnd_int(&seed, 0, 9));
            double w = X * (0.05 + rnd_unit(&seed) * 0.6);
            if (rnd64(&seed) & 1) { X = floor(X) + 1; w = floor(w) + 1; }
            double h = w / sqrt(2.0);
            double ox = (rnd64(&seed) & 1) ? 0 : (double)rnd_int(&seed, -500, 500);
            double oy = (rnd64(&seed) & 1) ? 0 : (double)rnd_int(&seed, -500, 500);
            c[0].x = ox - X; c[0].y = oy;
            c[1].x = ox - X + w; c[1].y = oy + h;
            c[2].x = ox + X - w; c[2].y = oy + h;
            c[3].x = ox + X; c[3].y = oy;
        } else if (!strcmp(family, "ess")) {
            /* point-symmetric about q: c0 = q - 3a, c3 = q + 3a?? the tangents
             * are then parallel; instead: e2 = -e1 needs q - c0 and q - c3
             * opposite: take c0 = q - 3a, c3 = q + 3b with b nearly a */
            V q = { (double)rnd_int(&seed, -300, 300), (double)rnd_int(&seed, -300, 300) };
            V a = { (double)rnd_int(&seed, 5, 60), (double)rnd_int(&seed, -60, 60) };
            V b = { a.x, a.y + (double)rnd_int(&seed, 1, 3) };
            int j = (int)rnd_int(&seed, 0, 1) * 3;
            c[0] = vsub(q, vmul(a, 3.0));
            c[3] = vadd(q, vmul(b, 3.0));
            c[1] = vadd(c[0], vmul(a, (double)j));
            c[2] = vsub(c[3], vmul(b, (double)j));
        } else if (!strcmp(family, "near")) {
            long size = (long[]){ 30, 100, 400, 1000 }[rnd64(&seed) % 4];
            V p0 = { (double)rnd_int(&seed, -size, size), (double)rnd_int(&seed, -size, size) };
            V q = { (double)rnd_int(&seed, -size, size), (double)rnd_int(&seed, -size, size) };
            V p3 = { (double)rnd_int(&seed, -size, size), (double)rnd_int(&seed, -size, size) };
            long j = (long[]){ 1, 3, 10, 30 }[rnd64(&seed) % 4];
            c[0] = p0; c[3] = p3;
            c[1].x = p0.x + rint(2 * (q.x - p0.x) / 3) + (double)rnd_int(&seed, -j, j);
            c[1].y = p0.y + rint(2 * (q.y - p0.y) / 3) + (double)rnd_int(&seed, -j, j);
            c[2].x = p3.x + rint(2 * (q.x - p3.x) / 3) + (double)rnd_int(&seed, -j, j);
            c[2].y = p3.y + rint(2 * (q.y - p3.y) / 3) + (double)rnd_int(&seed, -j, j);
        } else {
            fprintf(stderr, "unknown family %s\n", family);
            return 1;
        }
        V e1, e2, q1;
        if (!single_errors(c, 1e300, &e1, &e2, &q1)) continue;
        long double M = max_error_l(e1, e2, 0);
        if (!(M > 1e-9L) || !(M < 1e12L)) continue;
        double T = (double)M;
        for (int k = 0; k < 3; k++) T = nextafter(T, 0.0);
        for (int k = -3; k <= 3; k++, T = nextafter(T, INFINITY)) emit(c, T);
        made++;
    }
    fprintf(stderr, "%s: %ld cubics, %ld lines\n", family, made, made * 7);
    return 0;
}
