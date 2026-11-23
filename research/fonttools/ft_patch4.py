#!/usr/bin/env python3
"""Fourth change to fontTools/cu2qu/cu2qu.py: the fit test stops after 64
halvings.

Usage: ft_patch4.py <path to cu2qu.py>
"""
import sys

p = sys.argv[1]
s = open(p).read()


def rep(a, b):
    global s
    assert s.count(a) == 1, (s.count(a), a[:70])
    s = s.replace(a, b)


rep('''@cython.locals(mid=cython.complex, deriv3=cython.complex)
def cubic_farthest_fit_inside(p0, p1, p2, p3, tolerance):
    """Check if a cubic Bezier lies within a given distance of the origin.

    "Origin" means *the* origin (0,0), not the start of the curve. Note that no
    checks are made on the start and end positions of the curve; this function
    only checks the inside of the curve.

    Args:
        p0 (complex): Start point of curve.
        p1 (complex): First handle of curve.
        p2 (complex): Second handle of curve.
        p3 (complex): End point of curve.
        tolerance (double): Distance from origin.

    Returns:
        bool: True if the cubic Bezier ``p`` entirely lies within a distance
        ``tolerance`` of the origin, False otherwise.
    """
    # First check p2 then p1, as p2 has higher error early on.
    if abs(p2) <= tolerance and abs(p1) <= tolerance:
        return True

    # Split. The test is written for a NaN to end the recursion as well.
    mid = (p0 + 3 * (p1 + p2) + p3) * 0.125
    if not abs(mid) <= tolerance:
        return False
    deriv3 = (p3 + p2 - p1 - p0) * 0.125
    return cubic_farthest_fit_inside(
        p0, (p0 + p1) * 0.5, mid - deriv3, mid, tolerance
    ) and cubic_farthest_fit_inside(mid, mid + deriv3, (p2 + p3) * 0.5, p3, tolerance)
''', '''@cython.locals(mid=cython.complex, deriv3=cython.complex, depth=cython.int)
def cubic_farthest_fit_inside(p0, p1, p2, p3, tolerance, depth):
    """Check if a cubic Bezier lies within a given distance of the origin.

    "Origin" means *the* origin (0,0), not the start of the curve. Note that no
    checks are made on the start and end positions of the curve; this function
    only checks the inside of the curve.

    Args:
        p0 (complex): Start point of curve.
        p1 (complex): First handle of curve.
        p2 (complex): Second handle of curve.
        p3 (complex): End point of curve.
        tolerance (double): Distance from origin.
        depth (int): How many times the curve tested was halved to get this
            piece of it: 0 for the curve itself.

    Returns:
        bool: True if the cubic Bezier ``p`` entirely lies within a distance
        ``tolerance`` of the origin, False otherwise.
    """
    # First check p2 then p1, as p2 has higher error early on.
    if abs(p2) <= tolerance and abs(p1) <= tolerance:
        return True

    # Split. The test is written for a NaN to end the recursion as well.
    mid = (p0 + 3 * (p1 + p2) + p3) * 0.125
    if not abs(mid) <= tolerance:
        return False
    # A curve that touches the tolerance can be split without end: a control
    # point stays the last bit outside while the middle stays inside, and one
    # half is the same curve again. After this many halvings the points of a
    # piece differ by rounding alone, and the piece is where its middle is.
    if depth == MAX_HALVINGS:
        return True
    deriv3 = (p3 + p2 - p1 - p0) * 0.125
    return cubic_farthest_fit_inside(
        p0, (p0 + p1) * 0.5, mid - deriv3, mid, tolerance, depth + 1
    ) and cubic_farthest_fit_inside(
        mid, mid + deriv3, (p2 + p3) * 0.5, p3, tolerance, depth + 1
    )
''')

rep('''    if not cubic_farthest_fit_inside(0, c1 - cubic[1], c2 - cubic[2], 0, tolerance):
        return None
''', '''    if not cubic_farthest_fit_inside(0, c1 - cubic[1], c2 - cubic[2], 0, tolerance, 0):
        return None
''')

rep('''MAX_N = 100
''', '''MAX_N = 100

# How many times the fit test halves a curve at most.
MAX_HALVINGS = 64
''')
open(p, "w").write(s)
