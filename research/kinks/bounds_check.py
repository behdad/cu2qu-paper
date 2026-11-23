"""How far from 2/3 of the way to the tangents' intersection can a handle be,
for the single quadratic to pass cu2qu's fit test?

The quadratic, raised to a cubic, has control points c1' = c0 + (2/3)(q - c0)
and c2' = c3 + (2/3)(q - c3); with e1 = c1' - c1 and e2 = c2' - c2 the two
curves differ, at equal parameters, by

    E(t) = 3 t (1-t) ((1-t) e1 + t e2)

and the test passes when max |E| <= tolerance.  e1 is along the first handle
and e2 along the second.  This computes, by brute force, the largest |e1| that
any e2 on a line at angle psi to e1 allows, as a multiple of the tolerance;
the claims are 2 sqrt(3) = 3.4641 for psi -> 0 (e2 = -e1) and 9/4 for psi =
90 degrees (e2 = 0), and nothing above 2 sqrt(3) anywhere.

Then it checks the limit against cu2qu itself: cubics built with a given
(e1, e2) are converted at tolerance 1 with |e1| just below and just above the
computed limit.
"""

import cmath
import math

import numpy as np

from fontTools.cu2qu.cu2qu import curve_to_quadratic

T = np.linspace(0.0, 1.0, 20001)
B1 = 3 * T * (1 - T) ** 2
B2 = 3 * T**2 * (1 - T)


def max_err(e1, e2):
    return float(np.max(np.abs(B1 * e1 + B2 * e2)))


def best(psi):
    """Smallest max |E| for e1 = 1 over e2 = lam * exp(i psi)."""
    u = cmath.rect(1.0, psi)
    lo, hi = -3.0, 3.0
    lams = np.linspace(lo, hi, 6001)
    vals = [max_err(1.0, lam * u) for lam in lams]
    i = int(np.argmin(vals))
    # refine
    a, b = lams[max(i - 1, 0)], lams[min(i + 1, len(lams) - 1)]
    for _ in range(60):
        m1 = a + (b - a) / 3
        m2 = b - (b - a) / 3
        if max_err(1.0, m1 * u) < max_err(1.0, m2 * u):
            b = m2
        else:
            a = m1
    lam = (a + b) / 2
    return max_err(1.0, lam * u), lam


def cubic_with(e1, e2, q0, q1, q2):
    """The cubic whose single quadratic is q0 q1 q2 with these errors."""
    c1 = q0 + (q1 - q0) * (2 / 3) - e1
    c2 = q2 + (q1 - q2) * (2 / 3) - e2
    return [(z.real, z.imag) for z in (q0, c1, c2, q2)]


def main():
    print("psi (degrees between the end tangents)   largest |e1| / tolerance   at e2 / e1")
    worst = 0.0
    for deg in (0.001, 1, 5, 10, 20, 30, 45, 60, 75, 90, 120, 150, 179.999):
        m, lam = best(math.radians(deg))
        c = 1 / m
        worst = max(worst, c)
        print("%10.3f %32.4f %16.4f" % (deg, c, lam))
    print("largest over all: %.4f; 2 sqrt(3) = %.4f; 9/4 = %.4f" % (worst, 2 * math.sqrt(3), 2.25))

    print()
    print("against cu2qu, tolerance 1: quadratics for |e1| at 0.999 and 1.001 of the limit")
    for deg in (10, 45, 90, 150):
        psi = math.radians(deg)
        m, lam = best(psi)
        c = 1 / m
        # A quadratic whose end tangents are psi apart, with handles of 40.
        q0 = 0j
        q1 = 40 + 0j
        d2 = cmath.rect(1.0, psi)  # direction of e2: along q1 - q2
        q2 = q1 - 40 * d2
        out = []
        for f in (0.999, 1.001):
            e1 = f * c * (q1 - q0) / abs(q1 - q0)
            e2 = f * c * lam * d2
            cubic = cubic_with(e1, e2, q0, q1, q2)
            out.append(len(curve_to_quadratic(cubic, 1.0)) - 2)
        print("  psi %5.1f  limit %.4f  quadratics %s" % (deg, c, out))


if __name__ == "__main__":
    main()
