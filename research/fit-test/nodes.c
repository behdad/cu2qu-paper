/* How much work the recursive fit test does at its limit, with the depth
 * limit that kurbo has (64 halvings, then true): for families of error
 * curves (0, e1, e2, 0) and tolerances within K doubles of the largest
 * distance, count the pieces tested, and how deep the test went.
 *
 * Usage: nodes <family> <seed> <count> [K]
 *   generic   e1, e2 random, of comparable size
 *   line      e2 = 0 (the error curve is on a line; farthest at t = 1/3)
 *   small     e2 tiny against e1 (nearly on a line)
 *   axis      e1 = s (p, 1), e2 = s (0, -2): farthest on the x axis, at t = 1/3
 *   flat      |e1| = |e2|, (e2 - e1) at right angles to (e1 + e2) and
 *             sqrt(2) times as long: the distance is flat to fourth order
 *             at the middle
 *   flatrot   the same, rotated by a random angle
 *   nearflat  flat, with the length ratio off by up to 1e-3
 *   ess       e2 = -e1 or nearly: two places equally far
 * Prints the distribution of pieces tested and of the depth reached.
 */
#include <stdio.h>
#include <stdlib.h>
#include "fitsim.h"

int main(int argc, char **argv)
{
    const char *family = argv[1];
    uint64_t seed = strtoull(argv[2], 0, 10) * 0x9E3779B97F4A7C15ULL + 777;
    long count = atol(argv[3]);
    int K = argc > 4 ? atoi(argv[4]) : 3;
    for (int i = 0; i < 20; i++) rnd64(&seed);
    long hist_nodes[64], hist_depth[80], capped = 0, tests = 0, res[2] = { 0, 0 };
    memset(hist_nodes, 0, sizeof hist_nodes);
    memset(hist_depth, 0, sizeof hist_depth);
    long worst = 0;
    V w1 = { 0, 0 }, w2 = { 0, 0 };
    double wt = 0;
    for (long n = 0; n < count; n++) {
        V e1, e2;
        double s = ldexp(1.0 + rnd_unit(&seed), (int)rnd_int(&seed, -3, 6));
        double ang = rnd_unit(&seed) * 6.283185307179586;
        if (!strcmp(family, "generic")) {
            e1.x = (rnd_unit(&seed) * 2 - 1) * s; e1.y = (rnd_unit(&seed) * 2 - 1) * s;
            e2.x = (rnd_unit(&seed) * 2 - 1) * s; e2.y = (rnd_unit(&seed) * 2 - 1) * s;
        } else if (!strcmp(family, "line")) {
            e1.x = s * cos(ang); e1.y = (rnd64(&seed) & 1) ? 0.0 : s * sin(ang);
            e2.x = 0; e2.y = 0;
            if (rnd64(&seed) & 1) { V t = e1; e1 = e2; e2 = t; }
        } else if (!strcmp(family, "small")) {
            e1.x = s * cos(ang); e1.y = s * sin(ang);
            double tiny = s * ldexp(1.0, -(int)rnd_int(&seed, 20, 60));
            e2.x = (rnd_unit(&seed) * 2 - 1) * tiny; e2.y = (rnd_unit(&seed) * 2 - 1) * tiny;
            if (rnd64(&seed) & 1) { V t = e1; e1 = e2; e2 = t; }
        } else if (!strcmp(family, "axis")) {
            double p = 1.5 + rnd_unit(&seed) * 30;
            e1.x = s * p; e1.y = s; e2.x = 0; e2.y = -2 * s;
            if (rnd64(&seed) & 1) { double t = e1.x; e1.x = e1.y; e1.y = t; t = e2.x; e2.x = e2.y; e2.y = t; }
        } else if (!strncmp(family, "flat", 4) || !strcmp(family, "nearflat")) {
            double ratio = sqrt(2.0);
            if (!strcmp(family, "nearflat")) ratio *= 1 + (rnd_unit(&seed) * 2 - 1) * ldexp(1.0, -(int)rnd_int(&seed, 10, 40));
            double c = 1, sn = 0;
            if (!strcmp(family, "flatrot") || (!strcmp(family, "nearflat") && (rnd64(&seed) & 1))) { c = cos(ang); sn = sin(ang); }
            /* m = s (0, 1), r = s (ratio, 0), rotated */
            V m = { -sn * s, c * s }, r = { c * s * ratio, sn * s * ratio };
            e1 = vsub(m, r); e2 = vadd(m, r);
        } else if (!strcmp(family, "ess")) {
            e1.x = s * cos(ang); e1.y = s * sin(ang);
            e2.x = -e1.x; e2.y = -e1.y;
            if (rnd64(&seed) & 1) { e2.x *= 1 + ldexp(rnd_unit(&seed), -20); e2.y *= 1 - ldexp(rnd_unit(&seed), -20); }
        } else {
            fprintf(stderr, "unknown family\n");
            return 1;
        }
        long double M = max_error_l(e1, e2, 0);
        double T = (double)M;
        for (int k = 0; k < K; k++) T = nextafter(T, 0.0);
        for (int k = -K; k <= K; k++, T = nextafter(T, INFINITY)) {
            int r = fit_cap(e1, e2, T, 64, 200000000L);
            tests++;
            if (r == FIT_DEEP) { capped++; continue; }
            res[r]++;
            int b = 0;
            while ((1L << b) < fit_nodes) b++;
            hist_nodes[b]++;
            hist_depth[fit_depth]++;
            if (fit_nodes > worst) { worst = fit_nodes; w1 = e1; w2 = e2; wt = T; }
        }
    }
    printf("%s: %ld tests (tolerance within %d doubles of the largest distance): true %ld, false %ld, gave up after 2e8 pieces %ld\n",
           family, tests, K, res[1], res[0], capped);
    printf("  pieces tested, by power of two:");
    for (int b = 0; b < 64; b++) if (hist_nodes[b]) printf(" <=%ld:%ld", 1L << b, hist_nodes[b]);
    printf("\n  depth reached:");
    for (int d = 0; d < 80; d++) if (hist_depth[d]) printf(" %d:%ld", d, hist_depth[d]);
    printf("\n  most pieces: %ld for e1 (%a, %a) e2 (%a, %a) tolerance %a\n", worst, w1.x, w1.y, w2.x, w2.y, wt);
    return 0;
}
