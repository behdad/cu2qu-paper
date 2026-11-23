#!/usr/bin/env python3
"""Third change to fontTools/cu2qu/cu2qu.py, the small form: the single
quadratic is not tried when the cubic's third difference rules it out. The
check lives in cubic_approx_quadratic, which both APIs go through.

Usage: ft_patch3b.py <path to cu2qu.py>
"""
import sys

p = sys.argv[1]
s = open(p).read()


def rep(a, b):
    global s
    assert s.count(a) == 1, (s.count(a), a[:70])
    s = s.replace(a, b)


rep('''@cython.cfunc
@cython.inline
@cython.locals(tolerance=cython.double)
@cython.locals(
    q1=cython.complex,
    c0=cython.complex,
    c1=cython.complex,
    c2=cython.complex,
    c3=cython.complex,
)
def cubic_approx_quadratic(cubic, tolerance):
    """Approximate a cubic Bezier with a single quadratic within a given tolerance.

    Args:
        cubic (sequence): Four complex numbers representing control points of
            the cubic Bezier curve.
        tolerance (double): Permitted deviation from the original curve.

    Returns:
        Three complex numbers representing control points of the quadratic
        curve if it fits within the given tolerance, or ``None`` if no suitable
        curve could be calculated.
    """

    q1 = calc_intersect(cubic[0], cubic[1], cubic[2], cubic[3])
''', '''@cython.cfunc
@cython.inline
@cython.returns(cython.double)
@cython.locals(
    p0=cython.complex, p1=cython.complex, p2=cython.complex, p3=cython.complex
)
def cubic_third_difference(cubic):
    """Return the length of the third difference of a cubic Bezier curve.

    Args:
        cubic (sequence): Four complex numbers representing control points of
            the cubic Bezier curve.

    Returns:
        double: ``abs(p3 - 3*p2 + 3*p1 - p0)``: zero for a cubic that is a
        quadratic, and what decides how many quadratics another one takes.
    """
    p0 = cubic[0]
    p1 = cubic[1]
    p2 = cubic[2]
    p3 = cubic[3]
    return abs(p3 - 3 * p2 + 3 * p1 - p0)


# A single quadratic is not tried for a third difference of more than this
# many times the tolerance. A quadratic with the end points of a cubic gets at
# least sqrt(3)/36 of the cubic's third difference away from it, wherever its
# off-curve point is, so none fits above 36/sqrt(3) = 20.78; this is past that
# by more than the fit test can be wrong.
QUADRATIC_MAX_D3 = 21.0


@cython.cfunc
@cython.inline
@cython.locals(tolerance=cython.double, limit=cython.double)
@cython.locals(
    q1=cython.complex,
    c0=cython.complex,
    c1=cython.complex,
    c2=cython.complex,
    c3=cython.complex,
)
def cubic_approx_quadratic(cubic, tolerance):
    """Approximate a cubic Bezier with a single quadratic within a given tolerance.

    Args:
        cubic (sequence): Four complex numbers representing control points of
            the cubic Bezier curve.
        tolerance (double): Permitted deviation from the original curve.

    Returns:
        Three complex numbers representing control points of the quadratic
        curve if it fits within the given tolerance, or ``None`` if no suitable
        curve could be calculated.
    """

    # Not worth trying if the third difference rules it out.
    limit = QUADRATIC_MAX_D3
    if not cubic_third_difference(cubic) <= tolerance * limit:
        return None

    q1 = calc_intersect(cubic[0], cubic[1], cubic[2], cubic[3])
''')

rep('''@cython.cfunc
@cython.returns(cython.int)
@cython.locals(tolerance=cython.double, d3=cython.double, limit=cython.double)
@cython.locals(n=cython.int)
@cython.locals(
    p0=cython.complex, p1=cython.complex, p2=cython.complex, p3=cython.complex
)
def cubic_spline_n(cubic, tolerance):''', '''@cython.cfunc
@cython.returns(cython.int)
@cython.locals(tolerance=cython.double, d3=cython.double, limit=cython.double)
@cython.locals(n=cython.int)
def cubic_spline_n(cubic, tolerance):''')

rep('''    p0 = cubic[0]
    p1 = cubic[1]
    p2 = cubic[2]
    p3 = cubic[3]
    d3 = abs(p3 - 3 * p2 + 3 * p1 - p0)
    for n in range(2, MAX_N + 1):
''', '''    d3 = cubic_third_difference(cubic)
    for n in range(2, MAX_N + 1):
''')
open(p, "w").write(s)
