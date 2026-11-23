"""Shared helpers for the research-error scripts."""
import math
import numpy as np

S = '/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad'
R = S + '/research-error'
W = '/home/behdad/hb/cu2qu-work/research-error'
NF = 8

def load_cubics(path=W + '/corpus.bin'):
    return np.fromfile(path, dtype=np.float64).reshape(-1, 4, 2)

def load_gt(path, nmax):
    a = np.fromfile(path, dtype=np.float64)
    return a.reshape(-1, nmax, NF)

def Phi_ramp(n):
    """Largest |phi| of the cu2qu ramp, n >= 2 (equal-parameter error is |d3| Phi / n^3)."""
    if n == 2:
        return 4 / 27.0
    g = (n - 2) / (4.0 * (n - 1))
    first = 4 / 27.0 * (1 + g) ** 3
    w2 = (0.25 - 0.5 / (n - 1)) / 3
    second = (n - 3) / (4.0 * (n - 1)) + 2 * w2 * math.sqrt(w2)
    return max(first, second)

def d3_of(c):
    return c[:, 3] - 3 * c[:, 2] + 3 * c[:, 1] - c[:, 0]

def first_fit(err, tol):
    """err: (N, nmax) array of errors for n = 1..nmax (NaN = not computed, inf = not available).
    Returns the least n with err <= tol, 0 if none."""
    ok = err <= tol
    n = np.argmax(ok, axis=1) + 1
    n[~ok.any(axis=1)] = 0
    return n
