/* List the font cubics whose single quadratic is nearest its limit.
 * Usage: neartie <cubics.txt> <tolerance> <relative threshold> */
#include <stdio.h>
#include <stdlib.h>
#include "fitsim.h"
#include "closed.h"

int main(int argc, char **argv)
{
    FILE *f = fopen(argv[1], "r");
    double tol = atof(argv[2]), thr = atof(argv[3]);
    double v[8];
    int i = 0;
    if (!f) return 1;
    for (;; i++) {
        int k;
        for (k = 0; k < 8; k++) if (fscanf(f, "%la", &v[k]) != 1) break;
        if (k < 8) break;
        V c[4];
        for (k = 0; k < 4; k++) { c[k].x = v[2 * k]; c[k].y = v[2 * k + 1]; }
        V e1, e2, q;
        if (!single_errors(c, tol, &e1, &e2, &q)) continue;
        long double M = max_error_l(e1, e2, 0);
        double rel = (double)(M / tol) - 1;
        if (fabs(rel) > thr) continue;
        int r = fit_cap(e1, e2, tol, 64, 1L << 40);
        printf("cubic %d: (%g,%g) (%g,%g) (%g,%g) (%g,%g)  e1 (%.17g,%.17g) e2 (%.17g,%.17g)  max/tol-1 = %.3g  recursive %d (%ld pieces, depth %d) closed %d\n",
               i, c[0].x, c[0].y, c[1].x, c[1].y, c[2].x, c[2].y, c[3].x, c[3].y, e1.x, e1.y, e2.x, e2.y, rel, r, fit_nodes, fit_depth, closed_fit(e1, e2, tol));
    }
    return 0;
}
