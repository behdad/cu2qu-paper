#!/usr/bin/env python3
"""What does fontTools' cu2qu do when a cubic's third difference is exactly
at a limit of the closed form (54 * tol for n = 2, 128 * tol for n = 3)?
Run with PYTHONPATH pointing at the fontTools to test."""
import sys
from fontTools.cu2qu import cu2qu
print('fontTools cu2qu from', cu2qu.__file__, 'compiled:', cu2qu.COMPILED)

def d3(c):
    z = [complex(*p) for p in c]
    return abs(z[3] - 3 * z[2] + 3 * z[1] - z[0])

cases = {
    'n=2 tie (|d3| = 54, tol 1)':  [(0, 0), (0, 100), (100, 200), (354, 300)],
    'n=2 tie, vertical d3':        [(0, 0), (100, 0), (200, 100), (300, 354)],
    'n=3 tie (|d3| = 128, tol 1)': [(0, 0), (0, 100), (100, 200), (428, 300)],
    'n=3 tie, vertical d3':        [(0, 0), (100, 0), (200, 100), (300, 428)],
    'just under n=2 limit':        [(0, 0), (0, 100), (100, 200), (353.999, 300)],
    'just over n=2 limit':         [(0, 0), (0, 100), (100, 200), (354.001, 300)],
}
for name, c in cases.items():
    try:
        sp = cu2qu.curve_to_quadratic(c, 1.0)
        res = f'n = {len(sp) - 2}'
    except RecursionError:
        res = 'RecursionError'
    except Exception as ex:
        res = type(ex).__name__
    print(f'{name:32s} |d3| = {d3(c):.6f}  ->  {res}')
