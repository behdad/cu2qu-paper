#!/usr/bin/env python3
"""Check the piece that cu2qu's fit test never gets past, in one coordinate:
(A, A + ulp, A, A) with a tolerance of A. Its middle is A (inside), its
second point is outside, and its first half is the same piece again,
when A's 53-bit significand is 3 mod 8 and A is at least 8/7 of a power of
two. Tested on random A with the arithmetic of split_cubic_into_two.

Usage: trap_check.py [count]
"""
import math, random, sys

count = int(sys.argv[1]) if len(sys.argv) > 1 else 2000000
rnd = random.Random(7)


def halves(p0, p1, p2, p3):
    mid = (p0 + 3 * (p1 + p2) + p3) * 0.125
    deriv3 = (p3 + p2 - p1 - p0) * 0.125
    return mid, (p0, (p0 + p1) * 0.5, mid - deriv3, mid), (mid, mid + deriv3, (p2 + p3) * 0.5, p3)


by_class = {}
for _ in range(count):
    m = rnd.getrandbits(52) | (1 << 52)
    A = math.ldexp(m, -52)  # in [1, 2)
    up = math.nextafter(A, math.inf)
    piece = (A, up, A, A)
    mid, left, right = halves(*piece)
    trapped = mid == A and left == piece
    predicted = m % 8 == 3 and A >= 8 / 7
    key = (m % 8, A >= 8 / 7)
    t = by_class.setdefault(key, [0, 0, 0])
    t[0] += 1
    t[1] += trapped
    t[2] += trapped != predicted
    # the mirrored piece
    piece2 = (A, A, up, A)
    mid2, left2, right2 = halves(*piece2)
    assert (mid2 == A and right2 == piece2) == trapped, (A, piece2)
print("significand mod 8, A >= 8/7: pieces, of which trapped, of which not as predicted")
for key in sorted(by_class):
    print(" ", key, by_class[key])
print("mispredicted:", sum(t[2] for t in by_class.values()), "of", count)
