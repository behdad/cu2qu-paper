/* What the fit test for the single quadratic does on cubics from fonts:
 * how often it is reached, how often its first level answers, how many
 * pieces it tests otherwise; and whether the closed form (closed.h) gives
 * the same answer.
 *
 * Usage: profile <cubics.txt> <tolerance>...
 * cubics.txt: one cubic per line, 8 coordinates as C hex floats.
 */
#include <stdio.h>
#include <stdlib.h>
#include "fitsim.h"
#include "closed.h"

int main(int argc, char **argv)
{
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 1; }
    static V cubics[200000][4];
    int n = 0;
    while (n < 200000) {
        double v[8];
        int k;
        for (k = 0; k < 8; k++) if (fscanf(f, "%la", &v[k]) != 1) break;
        if (k < 8) break;
        for (k = 0; k < 4; k++) { cubics[n][k].x = v[2 * k]; cubics[n][k].y = v[2 * k + 1]; }
        n++;
    }
    fclose(f);
    printf("%d cubics\n", n);
    for (int a = 2; a < argc; a++) {
        double tol = atof(argv[a]);
        long tried = 0, accept0 = 0, reject0 = 0, deeper = 0, deeper_true = 0, differ = 0;
        long nodes_total = 0, nodes_max = 0, depth_max = 0, iters_total = 0, iters_max = 0;
        long hist[40];
        memset(hist, 0, sizeof hist);
        double closest = 1e300;
        for (int i = 0; i < n; i++) {
            V e1, e2, q;
            if (!single_errors(cubics[i], tol, &e1, &e2, &q)) continue;
            tried++;
            int r = fit_cap(e1, e2, tol, 64, 1L << 40);
            int c = closed_fit(e1, e2, tol);
            if (r != c) {
                differ++;
                printf("  DIFFER cubic %d: recursive %d closed %d\n", i, r, c);
            }
            if (fit_nodes == 1) {
                if (r) accept0++; else reject0++;
            } else {
                deeper++;
                deeper_true += r;
                nodes_total += fit_nodes;
                if (fit_nodes > nodes_max) nodes_max = fit_nodes;
                if (fit_depth > depth_max) depth_max = fit_depth;
                iters_total += closed_iters;
                if (closed_iters > iters_max) iters_max = closed_iters;
                int b = 0;
                while ((1L << b) < fit_nodes) b++;
                hist[b]++;
                double rel = fabs((double)(max_error_l(e1, e2, 0) / tol) - 1);
                if (rel < closest) closest = rel;
            }
        }
        printf("tolerance %g: single quadratic tried for %ld of %d (%.1f%%)\n", tol, tried, n, 100.0 * tried / n);
        printf("  first level: fits %ld (%.1f%%), does not %ld (%.1f%%); goes deeper %ld (%.1f%%), of which %ld fit\n",
               accept0, 100.0 * accept0 / tried, reject0, 100.0 * reject0 / tried, deeper, 100.0 * deeper / tried, deeper_true);
        printf("  when deeper: %.1f pieces on average, at most %ld, depth at most %ld; pieces by power of two:",
               deeper ? (double)nodes_total / deeper : 0.0, nodes_max, depth_max);
        for (int b = 0; b < 40; b++) if (hist[b]) printf(" <=%ld:%ld", 1L << b, hist[b]);
        printf("\n  closed form on those: %.2f Newton steps on average, at most %ld; answers differing from the recursive test: %ld\n",
               deeper ? (double)iters_total / deeper : 0.0, iters_max, differ);
        printf("  closest that the largest distance of a deeper one comes to the tolerance: relative %.3g\n", closest);
    }
    return 0;
}
