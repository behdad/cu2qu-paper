"""Constructed masters that meet only the equal-ratio condition.

Two masters of one joint.  Coming in: a line of the same length L in both
masters.  Going out: a cubic whose first handle has length H in both masters
(so the ratio of the two handles at the joint is L/H in both), but whose
direction turns by theta between the masters, and whose handle reaches a
different fraction of the way to the intersection of the end tangents:
off-curve = on-curve + s * handle, with s = 1.5 (1 + x) in one master and
1.5 (1 - x) in the other.  The far end of each cubic is exactly as a
quadratic would have it, so the largest equal-parameter distance between the
cubic and its single quadratic is (4/9) x H.

The cubic masters interpolate without a kink (equal ratio).  The kink of the
interpolated quadratic masters is measured over the whole axis and compared
with

    tan(kink) = (R - 1) sin(theta) / (2 sqrt(R) + (1 + R) cos(theta))

with R the ratio of the two masters' s.

usage: adversarial.py [TOLERANCE]
"""

import cmath
import math
import sys

from fontTools.cu2qu.cu2qu import curves_to_quadratic

DEG = 180 / math.pi


def master(theta, H, s, L, turn=60.0, far=None):
    """Joint at the origin, tangent at angle theta (degrees)."""
    u = cmath.rect(1.0, theta / DEG)
    line = (-L * u, 0j)
    q1 = s * H * u  # where the end tangents meet
    v = cmath.rect(1.0, (theta + turn) / DEG)
    q2 = q1 + (far if far is not None else s * H) * v
    c1 = H * u
    c2 = q2 + (q1 - q2) * (2.0 / 3.0)
    return line, (0j, c1, c2, q2)


def angle(a, b):
    return cmath.phase(b * a.conjugate()) * DEG


def formula(R, theta):
    t = theta / DEG
    return math.atan2((R - 1) * math.sin(t), 2 * math.sqrt(R) + (1 + R) * math.cos(t)) * DEG


def run(theta, H, x, tol, L=30.0):
    sa, sb = 1.5 * (1 + x), 1.5 * (1 - x)
    la, ca = master(0.0, H, sa, L)
    lb, cb = master(theta, H, sb, L)
    as_tuples = lambda c: [(z.real, z.imag) for z in c]
    qa, qb = curves_to_quadratic([as_tuples(ca), as_tuples(cb)], [tol, tol])
    n = len(qa) - 2
    worst_q = worst_c = 0.0
    worst_dev = 0.0
    steps = 2000
    for i in range(1, steps):
        t = i / steps
        mix = lambda a, b: (1 - t) * a + t * b
        inh = mix(la[1] - la[0], lb[1] - lb[0])
        out_c = mix(ca[1] - ca[0], cb[1] - cb[0])
        out_q = mix(complex(*qa[1]) - ca[0], complex(*qb[1]) - cb[0])
        worst_c = max(worst_c, abs(angle(inh, out_c)))
        worst_q = max(worst_q, abs(angle(inh, out_q)))
        dev = abs((out_q * out_c.conjugate()).imag) / abs(out_c)
        worst_dev = max(worst_dev, dev)
    return n, worst_c, worst_q, abs(formula(sb / sa, theta)), worst_dev


def main():
    tol = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    print("tolerance %g; x is the largest that still gets a single quadratic: 2.25 tol / H" % tol)
    print(
        "%6s %7s %7s %3s  %12s %12s %12s %14s"
        % ("H", "theta", "x", "n", "cubic kink", "quad kink", "formula", "sideways dev")
    )
    for H in (5.0, 10.0, 20.0, 50.0, 100.0, 200.0):
        for theta in (5.0, 10.0, 20.0, 45.0, 90.0, 135.0):
            x = 2.2499 * tol / H
            if x >= 0.9:
                continue
            n, kc, kq, f, dev = run(theta, H, x, tol)
            print(
                "%6g %7g %7.4f %3d  %12.3e %12.4f %12.4f %14.4f"
                % (H, theta, x, n, kc, kq, f, dev)
            )
    print()
    print("just past the limit (x = 2.26 tol / H): no single quadratic, no kink")
    for H in (20.0, 100.0):
        n, kc, kq, f, dev = run(45.0, H, 2.26 * tol / H, tol)
        print("%6g %7g %3d  %12.3e %12.3e" % (H, 45.0, n, kc, kq))


if __name__ == "__main__":
    main()
