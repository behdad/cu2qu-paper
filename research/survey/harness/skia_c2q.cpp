// Runs Skia's own GrPathUtils::convertCubicToQuads (linked from a prebuilt
// libskia.a) over the corpus.
//
//   skia_c2q <corpus.bin> <tolScale> <out.bin> [time-reps]
//
// Corpus: "C2QCOR01", u32 count, count * 8 f64.
// Output: "C2QRES01", u32 count, per cubic u32 nseg, nseg * {u32 kind, 6 f64}.

#include "include/core/SkPoint.h"
#include "include/private/SkTArray.h"
#include "src/gpu/ganesh/geometry/GrPathUtils.h"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s corpus.bin tolScale out.bin [reps]\n", argv[0]);
        return 2;
    }
    FILE* f = fopen(argv[1], "rb");
    if (!f) { perror("corpus"); return 1; }
    char magic[8];
    uint32_t count;
    if (fread(magic, 1, 8, f) != 8 || memcmp(magic, "C2QCOR01", 8) != 0) return 1;
    if (fread(&count, 4, 1, f) != 1) return 1;
    std::vector<double> pts(size_t(count) * 8);
    if (fread(pts.data(), 8, pts.size(), f) != pts.size()) return 1;
    fclose(f);

    float tol = (float)atof(argv[2]);
    int reps = argc > 4 ? atoi(argv[4]) : 0;

    FILE* o = fopen(argv[3], "wb");
    if (!o) { perror("out"); return 1; }
    fwrite("C2QRES01", 1, 8, o);
    fwrite(&count, 4, 1, o);
    skia_private::TArray<SkPoint, true> quads;
    for (uint32_t i = 0; i < count; i++) {
        const double* c = &pts[size_t(i) * 8];
        SkPoint p[4] = {{(float)c[0], (float)c[1]}, {(float)c[2], (float)c[3]},
                        {(float)c[4], (float)c[5]}, {(float)c[6], (float)c[7]}};
        quads.clear();
        GrPathUtils::convertCubicToQuads(p, tol, &quads);
        uint32_t n = quads.size() / 3;
        fwrite(&n, 4, 1, o);
        for (uint32_t k = 0; k < n; k++) {
            uint32_t kind = 2;
            double v[6] = {quads[3 * k].fX, quads[3 * k].fY, quads[3 * k + 1].fX,
                           quads[3 * k + 1].fY, quads[3 * k + 2].fX, quads[3 * k + 2].fY};
            fwrite(&kind, 4, 1, o);
            fwrite(v, 8, 6, o);
        }
    }
    fclose(o);

    if (reps > 0) {
        double best = 1e30;
        size_t nq = 0;
        float check = 0;
        for (int r = 0; r < reps; r++) {
            auto t0 = std::chrono::steady_clock::now();
            nq = 0;
            for (uint32_t i = 0; i < count; i++) {
                const double* c = &pts[size_t(i) * 8];
                SkPoint p[4] = {{(float)c[0], (float)c[1]}, {(float)c[2], (float)c[3]},
                                {(float)c[4], (float)c[5]}, {(float)c[6], (float)c[7]}};
                quads.clear();
                GrPathUtils::convertCubicToQuads(p, tol, &quads);
                nq += quads.size() / 3;
                check += quads[quads.size() - 2].fX;
            }
            auto t1 = std::chrono::steady_clock::now();
            double dt = std::chrono::duration<double>(t1 - t0).count();
            if (dt < best) best = dt;
        }
        printf("{\"method\":\"skia_convertCubicToQuads\",\"tol\":%g,\"ns_per_cubic\":%.1f,\"quads\":%zu,\"check\":%g}\n",
               tol, best * 1e9 / count, nq, check);
    }
    return 0;
}
