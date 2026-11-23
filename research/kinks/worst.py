"""Show the joints where cu2qu changes the instance kink most.

usage: worst.py DATA.npz [COUNT] [exact|near|flag]
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stage2


def c(z):
    return "(%g,%g)" % (z.real, z.imag)


def main():
    path = sys.argv[1]
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    pop = sys.argv[3] if len(sys.argv) > 3 else "exact"
    data, meta, names = stage2.load(path)
    res = stage2.analyze(data, meta)
    mask = res["valid"] & {
        "exact": res["exact"],
        "near": res["near"],
        "flag": res["flag"] == 1,
    }[pop]
    idx = np.nonzero(mask)[0]
    order = idx[np.argsort(-res["dq"][idx])][:count]
    pts = data["pts"]
    kind = data["kind"]
    s1, s2, h1, h2 = stage2.single_fractions(pts)
    for k in order:
        j = res["j"][k]
        p = res["p"][k]
        print(
            "glyph %s  joint at segment %d: change by cu2qu %.3f deg (cubic instance %.3f, quadratic %.3f, rounded %.3f; rounding alone %.3f)  tangent turns %.2f deg between masters, at %s"
            % (
                names[data["gid"][j]],
                j,
                res["dq"][k],
                res["kc"][k],
                res["kq"][k],
                res["kr"][k],
                res["dr"][k],
                res["theta"][k],
                res["samples"][res["dq_loc"][k]],
            )
        )
        for m in range(pts.shape[1]):
            def seg(i):
                if kind[i] == 0:
                    return "line %s %s" % (c(pts[i, m, 0]), c(pts[i, m, 1]))
                return "cubic %s n=%d s=(%.4f,%.4f) |h|=(%.1f,%.1f)" % (
                    " ".join(c(z) for z in pts[i, m]),
                    data["n"][i],
                    s1[i, m],
                    s2[i, m],
                    abs(h1[i, m]),
                    abs(h2[i, m]),
                )
            print("   master %d  in: %s" % (m, seg(p)))
            print("             out: %s" % seg(j))
        print()


if __name__ == "__main__":
    main()
