"""Numerical check of the n >= 2 claims on real masters, with cu2qu's own code.

usage: check_linear.py DATA.npz [MAX_CURVES]

For every cubic that curves_to_quadratic gave n >= 2 quadratics (the same n in
every master):

 1. the spline of an interpolated cubic is the interpolation of the splines:
    cubic_spline(sum w_i C_i, n) == sum w_i cubic_spline(C_i, n), for weights
    that sum to one, some of them negative;
 2. the first and last off-curve points are 3/(2n) of the cubic's handles from
    the ends;
 3. the instance needs no more quadratics than the masters got, for weights
    that are not negative (cubic_spline_n of the interpolated cubic <= n).

Differences are reported relative to the size of the curve (the largest
coordinate difference among the control points of the masters), in units of
2**-52.
"""

import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stage2

from fontTools.cu2qu import cu2qu

assert not cu2qu.COMPILED, "needs the pure Python module, to call its internals"

EPS = 2.0**-52


def main():
    path = sys.argv[1]
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 10**9
    data, meta, names = stage2.load(path)
    pts = data["pts"]
    n = data["n"]
    tol = meta["tol"]
    S, M, _ = pts.shape
    rng = random.Random(1)
    idx = np.nonzero((data["kind"] == 1) & (n >= 2))[0]
    if len(idx) > limit:
        idx = np.array(sorted(rng.sample(list(idx), limit)))
    worst_lin = worst_handle = 0.0
    more = 0
    checked = 0
    for s in idx:
        cubics = [[complex(z) for z in pts[s, m]] for m in range(M)]
        k = int(n[s])
        splines = [cu2qu.cubic_spline(c, k) for c in cubics]
        flat = [z for c in cubics for z in c]
        size = max(
            max(z.real for z in flat) - min(z.real for z in flat),
            max(z.imag for z in flat) - min(z.imag for z in flat),
            1e-300,
        )
        mag = max(max(abs(z.real), abs(z.imag)) for z in flat)
        for m in range(M):
            c = cubics[m]
            a = splines[m][1] - c[0] - (c[1] - c[0]) * (1.5 / k)
            b = c[3] - splines[m][-2] - (c[3] - c[2]) * (1.5 / k)
            worst_handle = max(worst_handle, abs(a) / mag, abs(b) / mag)
        for trial in range(3):
            if trial == 0:
                w = [rng.random() for _ in range(M)]
            else:
                w = [rng.uniform(-1, 2) for _ in range(M)]
            total = sum(w)
            if abs(total) < 0.1:
                continue
            w = [x / total for x in w]
            mix = [sum(w[m] * cubics[m][i] for m in range(M)) for i in range(4)]
            direct = cu2qu.cubic_spline(mix, k)
            interp = [
                sum(w[m] * splines[m][i] for m in range(M)) for i in range(k + 2)
            ]
            scale = mag * max(1.0, sum(abs(x) for x in w))
            d = max(abs(x - y) for x, y in zip(direct, interp)) / scale
            worst_lin = max(worst_lin, d)
            if trial == 0:
                need = cu2qu.cubic_spline_n(mix, tol)
                if need > k or need == 0:
                    more += 1
        checked += 1
    print(
        "%s: %d cubics with n >= 2 checked, %d masters" % (meta["name"], checked, M)
    )
    print(
        "  spline of interpolated cubic vs interpolated splines: largest difference %.1f ulp of the coordinates"
        % (worst_lin / EPS)
    )
    print(
        "  end handles vs 3/(2n) of the cubic's: largest difference %.1f ulp of the coordinates"
        % (worst_handle / EPS)
    )
    print("  instances (weights >= 0) that take more quadratics than the masters got: %d" % more)


if __name__ == "__main__":
    main()
