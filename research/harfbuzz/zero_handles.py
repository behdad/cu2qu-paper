#!/usr/bin/env python3
"""Zero-length handles in real outlines: where are they, and what does
cu2qu's spline do to the tangent there?

For every cubic with a control point on its end point, in the first glyphs
of some CFF fonts:
  - the cubic's own tangent at that end is toward the other control point;
  - the neighbouring segment's tangent at the shared point tells whether the
    outline is smooth there (the designer meant a smooth join) or a corner;
  - cu2qu's spline starts along the direction from the end point to its
    second off-curve point (its first off-curve point is on the end point),
    or, for a single quadratic, toward the quadratic's control point.

Usage: PYTHONPATH=<fontTools Lib> zero_handles.py max_glyphs tolerance font...
"""
import math, sys
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import RecordingPen
from fontTools.cu2qu.cu2qu import curve_to_quadratic

max_glyphs = int(sys.argv[1])
tol = float(sys.argv[2])


def angle(u, v):
    """Unsigned angle between two directions, in degrees."""
    cr = u[0] * v[1] - u[1] * v[0]
    dt = u[0] * v[0] + u[1] * v[1]
    return abs(math.degrees(math.atan2(cr, dt)))


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def start_dir(seg):
    """Tangent direction at the start of a segment (a list of points)."""
    for p in seg[1:]:
        d = sub(p, seg[0])
        if d != (0, 0):
            return d
    return None


def end_dir(seg):
    for p in reversed(seg[:-1]):
        d = sub(seg[-1], p)
        if d != (0, 0):
            return d
    return None


total_cubics = 0
ends = []   # (join angle with neighbour, spline direction error, n, size)
for path in sys.argv[3:]:
    font = TTFont(path, fontNumber=0)
    gs = font.getGlyphSet()
    for name in font.getGlyphOrder()[:max_glyphs]:
        pen = RecordingPen()
        try:
            gs[name].draw(pen)
        except Exception:
            continue
        contour = []
        cur = start = None
        contours = []
        for op, args in pen.value:
            if op == "moveTo":
                cur = start = args[0]
                contour = []
            elif op == "lineTo":
                if args[0] != cur:
                    contour.append([cur, args[0]])
                cur = args[0]
            elif op == "curveTo":
                if len(args) == 3:
                    contour.append([cur, args[0], args[1], args[2]])
                cur = args[-1]
            elif op in ("closePath", "endPath"):
                if op == "closePath" and cur != start and start is not None:
                    contour.append([cur, start])
                if contour:
                    contours.append((contour, op == "closePath"))
                contour = []
        for contour, closed in contours:
            n = len(contour)
            for i, seg in enumerate(contour):
                if len(seg) != 4:
                    continue
                total_cubics += 1
                c0, c1, c2, c3 = seg
                for which in ("start", "end"):
                    if which == "start" and c1 != c0:
                        continue
                    if which == "end" and c2 != c3:
                        continue
                    if c0 == c1 and c2 == c3:
                        continue          # a straight line drawn as a cubic
                    # the cubic's own tangent at this end
                    own = start_dir(seg) if which == "start" else end_dir(seg)
                    # the neighbour's tangent at the shared point
                    j = i - 1 if which == "start" else i + 1
                    if not closed and (j < 0 or j >= n):
                        continue
                    nb = contour[j % n]
                    other = end_dir(nb) if which == "start" else start_dir(nb)
                    if own is None or other is None:
                        continue
                    join = angle(other, own)
                    # what the spline does
                    spline = curve_to_quadratic(seg, tol)
                    k = len(spline) - 2
                    if which == "start":
                        pts = spline
                    else:
                        pts = spline[::-1]
                    # direction in which the spline leaves this end
                    leave = None
                    for p in pts[1:]:
                        d = sub(p, pts[0])
                        if d != (0, 0):
                            leave = d
                            break
                    want = own if which == "start" else (-own[0], -own[1])
                    err = angle(leave, want)
                    size = math.hypot(*sub(c3, c0))
                    ends.append((join, err, k, size))

print(f"cubics {total_cubics}; ends with a zero-length handle (not straight lines): {len(ends)}")
smooth = [e for e in ends if e[0] < 1.0]
nearly = [e for e in ends if 1.0 <= e[0] < 5.0]
corner = [e for e in ends if e[0] >= 5.0]
print(f"  at a smooth join (neighbour within 1 degree of the cubic's own tangent): {len(smooth)}")
print(f"  nearly smooth (1 to 5 degrees): {len(nearly)}")
print(f"  at a corner (5 degrees or more): {len(corner)}")


def stats(name, es):
    if not es:
        return
    errs = sorted(e[1] for e in es)
    q = lambda f: errs[min(len(errs) - 1, int(f * len(errs)))]
    over1 = sum(1 for x in errs if x > 1)
    over5 = sum(1 for x in errs if x > 5)
    over10 = sum(1 for x in errs if x > 10)
    ns = sorted(e[2] for e in es)
    print(f"  {name}: spline leaves off the cubic's tangent by median {q(0.5):.1f}, 90% {q(0.9):.1f}, max {errs[-1]:.1f} degrees; over 1: {over1}, over 5: {over5}, over 10: {over10}; quadratics median {ns[len(ns) // 2]}")


stats("smooth joins", smooth)
stats("nearly smooth", nearly)
stats("corners", corner)
