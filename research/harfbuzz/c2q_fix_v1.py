#!/usr/bin/env python3
# Candidate fix V1 for src/hb-gpu-cu2qu.hh, applied to the pristine file:
# a zero-length end tangent takes its direction from the other control point.
import sys
p = sys.argv[1]
s = open(p).read()
a = '''  double ax = c1.x - c0.x, ay = c1.y - c0.y;
  double dx = c3.x - c2.x, dy = c3.y - c2.y;

'''
b = '''  double ax = c1.x - c0.x, ay = c1.y - c0.y;
  double dx = c3.x - c2.x, dy = c3.y - c2.y;

  /* A control point on its end point leaves that end no tangent to
   * intersect.  The curve leaves the end toward the other control point,
   * so take the direction from that one: the quadratic's control point
   * is then that control point. */
  if (ax == 0 && ay == 0)
  {
    ax = c2.x - c0.x;
    ay = c2.y - c0.y;
  }
  if (dx == 0 && dy == 0)
  {
    dx = c3.x - c1.x;
    dy = c3.y - c1.y;
  }

'''
assert s.count(a) == 1
s = s.replace(a, b)
open(p, 'w').write(s)
