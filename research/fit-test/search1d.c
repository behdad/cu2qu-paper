/* Search for finite inputs on which fontTools' fit test for the single
 * quadratic recurses without end, among cubics with three equal points
 * (the tangents do not intersect; cu2qu takes the common point as the
 * quadratic's off-curve point): the error curve is then (0, e, 0, 0) or
 * (0, 0, e, 0), on one line, and its largest distance is 4/9 |e|.
 *
 * Usage: search1d <seed> <count> <mode> [K]
 *   mode 0: (X,0) (0,0) (0,0) (0,0)      error on the x axis: one dimension
 *   mode 1: (0,0) (0,0) (0,0) (X,0)
 *   mode 2: (X,Y) (0,0) (0,0) (0,0)      error on a slanted line
 *   mode 3: (0,0) (0,0) (0,0) (X,Y)
 *   mode 4: integers: (X,0) 0 0 0 with X an integer up to 4096 times a power of 2
 *   mode 5: (X,Y) integers
 * The tolerance is the double nearest to the largest distance, and the K
 * doubles on either side of it (default 4).
 * Prints one line per input that passes 300 halvings.
 */
#include <stdio.h>
#include <stdlib.h>
#include "fitsim.h"

int main(int argc, char **argv)
{
    uint64_t seed = strtoull(argv[1], 0, 10) * 0x9E3779B97F4A7C15ULL + 12345;
    long count = atol(argv[2]);
    int mode = atoi(argv[3]);
    int K = argc > 4 ? atoi(argv[4]) : 4;
    long found = 0, tried = 0, results[3] = { 0, 0, 0 };
    long depth_hist[400];
    memset(depth_hist, 0, sizeof depth_hist);
    for (int i = 0; i < 20; i++) rnd64(&seed);
    for (long n = 0; n < count; n++) {
        V c[4] = { { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 } };
        double X, Y = 0;
        if (mode == 4 || mode == 5) {
            X = (double)rnd_int(&seed, 1, 4096);
            if (mode == 5) Y = (double)rnd_int(&seed, -4096, 4096);
            if (rnd64(&seed) & 1) { X /= 64; Y /= 64; }
        } else {
            X = 1.0 + rnd_unit(&seed);
            if (mode == 2 || mode == 3) Y = (rnd_unit(&seed) * 2 - 1) * 2;
        }
        int which = (mode == 0 || mode == 2) ? 0 : (mode == 1 || mode == 3) ? 3 : (int)(rnd64(&seed) & 1) * 3;
        c[which].x = X; c[which].y = Y;
        V e1, e2, q;
        /* the tolerance is not known yet: use a huge one for the pre-checks */
        if (!single_errors(c, 1e300, &e1, &e2, &q)) continue;
        long double M = max_error_l(e1, e2, 0);
        double T0 = (double)M;
        double T = T0;
        for (int k = 0; k < K; k++) T = nextafter(T, 0.0);
        for (int k = -K; k <= K; k++, T = nextafter(T, INFINITY)) {
            V f1, f2;
            if (!single_errors(c, T, &f1, &f2, &q)) continue;
            int r = fit(f1, f2, T, 300);
            tried++;
            results[r + 1]++;
            if (r != FIT_DEEP) depth_hist[fit_depth]++;
            if (r == FIT_DEEP) {
                found++;
                if (found <= 40)
                    printf("DEEP mode %d k %+d tol %a (%.17g) cubic (%.17g,%.17g) (%.17g,%.17g) (%.17g,%.17g) (%.17g,%.17g) e1 (%a,%a) e2 (%a,%a) trap (%a,%a) (%a,%a) (%a,%a) (%a,%a)\n",
                           mode, k, T, T, c[0].x, c[0].y, c[1].x, c[1].y, c[2].x, c[2].y, c[3].x, c[3].y,
                           f1.x, f1.y, f2.x, f2.y,
                           fit_trap[0].x, fit_trap[0].y, fit_trap[1].x, fit_trap[1].y,
                           fit_trap[2].x, fit_trap[2].y, fit_trap[3].x, fit_trap[3].y);
            }
        }
    }
    int maxd = 0;
    for (int i = 0; i < 400; i++) if (depth_hist[i]) maxd = i;
    printf("mode %d: %ld tests, %ld without end, true %ld false %ld, deepest that ended %d\n",
           mode, tried, results[0], results[2], results[1], maxd);
    return 0;
}
