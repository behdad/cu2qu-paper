#!/usr/bin/env python3
"""Check a closed form for the error of the cu2qu spline.

Claim.  Cut a cubic into n >= 2 equal pieces and build the cu2qu spline.
With d3 = c3 - 3 c2 + 3 c1 - c0 the cubic's third difference and
e = d3 / n^3, the difference between quadratic k and its piece at equal
parameters u is

    E_k(u) = e * phi_k(u),

a multiple of the one vector e, where phi_k is a scalar cubic polynomial
that depends on n and k only:

    phi_k(u) = g_k - u (2 g_k + k/(n-1)) + u^2 (1 + g_k + g_{k+1} + k/(n-1)) - u^3
    g_k = (n - 2k) / (4 (n - 1))  for 0 < k < n,   g_0 = g_n = 0

(g_k e is the offset of the implied on-curve point k from the cubic).

So n pieces fit iff |d3| * Phi(n) <= tolerance * n^3, with
Phi(n) = max over k and u of |phi_k(u)|.
"""
import math, random, sys

def g(n, k):
    return 0.0 if k == 0 or k == n else (n - 2 * k) / (4.0 * (n - 1))

def phi_coeffs(n, k):
    kap = k / (n - 1.0)
    return (g(n, k), -(2 * g(n, k) + kap), 1 + g(n, k) + g(n, k + 1) + kap, -1.0)

def poly_max_abs(c):
    """max |c0 + c1 u + c2 u^2 + c3 u^3| on [0, 1]."""
    f = lambda u: c[0] + u * (c[1] + u * (c[2] + u * c[3]))
    cand = [0.0, 1.0]
    # derivative: c1 + 2 c2 u + 3 c3 u^2
    a, b, cc = 3 * c[3], 2 * c[2], c[1]
    disc = b * b - 4 * a * cc
    if abs(a) > 1e-300 and disc >= 0:
        for s in (1, -1):
            u = (-b + s * math.sqrt(disc)) / (2 * a)
            if 0 < u < 1:
                cand.append(u)
    return max(abs(f(u)) for u in cand)

def Phi(n):
    return max(poly_max_abs(phi_coeffs(n, k)) for k in range(n))

def Phi_closed(n):
    """The closed form: the first piece, or the piece next to it."""
    if n == 2:
        return 4 / 27.0
    g1 = (n - 2) / (4.0 * (n - 1))
    first = (4 / 27.0) * (1 + g1) ** 3
    eps = 1 / (2.0 * (n - 1))
    m1 = (n - 3) / (4.0 * (n - 1))
    w = math.sqrt(max(0.0, (0.25 - eps) / 3))
    second = max(m1 + 2 * w ** 3, g1)
    return max(first, second)

def n_from_formula(d3, tol, max_n=100):
    for n in range(2, max_n + 1):
        if d3 * Phi(n) <= tol * n ** 3:
            return n
    return None

# --- 1. The table, and the closed form against the direct maximum.
print(' n   Phi(n)     closed     n^3/Phi(n)')
for n in list(range(2, 13)) + [16, 24, 32, 48, 64, 100, 1000]:
    print(f'{n:4d} {Phi(n):.6f}  {Phi_closed(n):.6f}  {n**3 / Phi(n):14.2f}')
assert all(abs(Phi(n) - Phi_closed(n)) < 1e-12 for n in range(2, 2000))
assert all(Phi(n + 1) > Phi(n) for n in range(2, 2000))
print('closed form matches for n = 2..1999; Phi is increasing')
print('limit: 1/4 + sqrt(3)/36 =', 0.25 + math.sqrt(3) / 36, '; Phi(2) = 4/27 =', 4 / 27.0)

# --- 2. The identity itself, on random cubics, built the way cu2qu does.
def spline(c, n):
    """off-curve and on-curve points of the cu2qu spline, and the pieces."""
    c0, c1, c2, c3 = c
    b1 = 3 * (c1 - c0); b2 = 3 * (c0 - 2 * c1 + c2); b3 = c3 - 3 * c2 + 3 * c1 - c0
    B = lambda t: c0 + t * (b1 + t * (b2 + t * b3))
    D = lambda t: b1 + t * (2 * b2 + 3 * t * b3)
    pieces = []
    for k in range(n):
        t0, t1 = k / n, (k + 1) / n
        a, d = B(t0), B(t1)
        pieces.append((a, a + D(t0) / (3 * n), d - D(t1) / (3 * n), d))
    off = []
    for k, (a, b, cc, d) in enumerate(pieces):
        p1 = a + 1.5 * (b - a); p2 = d + 1.5 * (cc - d)
        off.append(p1 + (p2 - p1) * (k / (n - 1)))
    on = [c0] + [(off[k - 1] + off[k]) / 2 for k in range(1, n)] + [c3]
    return off, on, pieces

random.seed(1)
worst = 0.0
for trial in range(2000):
    c = [complex(random.uniform(-1000, 1000), random.uniform(-1000, 1000)) for _ in range(4)]
    n = random.randint(2, 40)
    d3 = c[3] - 3 * c[2] + 3 * c[1] - c[0]
    e = d3 / n ** 3
    off, on, pieces = spline(c, n)
    for k in range(n):
        a, b, cc, d = pieces[k]
        co = phi_coeffs(n, k)
        for j in range(0, 17):
            u = j / 16.0
            cubic = (1-u)**3 * a + 3*(1-u)**2*u * b + 3*(1-u)*u*u * cc + u**3 * d
            quad = (1-u)**2 * on[k] + 2*(1-u)*u * off[k] + u*u * on[k + 1]
            pred = e * (co[0] + u * (co[1] + u * (co[2] + u * co[3])))
            worst = max(worst, abs((quad - cubic) - pred) / max(1.0, abs(d3)))
print(f'identity E_k(u) = e phi_k(u): largest relative discrepancy over 2000 random cubics: {worst:.2e}')

# --- 3. Against fontTools: the n it finds by trial, when it is not 1.
try:
    from fontTools.cu2qu.cu2qu import curve_to_quadratic
except Exception as ex:
    print('fontTools not available:', ex); sys.exit(0)
random.seed(2)
same = diff = ones = 0
examples = []
for trial in range(200000):
    scale = 10 ** random.uniform(0.5, 3.5)
    c = [(random.uniform(-scale, scale), random.uniform(-scale, scale)) for _ in range(4)]
    tol = random.choice([0.25, 0.5, 1.0, 2.0])
    try:
        sp = curve_to_quadratic(c, tol)
    except Exception:
        continue
    n_ft = len(sp) - 2
    if n_ft == 1:
        ones += 1
        continue
    cc = [complex(*p) for p in c]
    d3 = abs(cc[3] - 3 * cc[2] + 3 * cc[1] - cc[0])
    n_f = n_from_formula(d3, tol)
    if n_f == n_ft: same += 1
    else:
        diff += 1
        if len(examples) < 5: examples.append((n_ft, n_f, d3 / tol, n_ft ** 3 / Phi(n_ft), (n_f or 0) and n_f ** 3 / Phi(n_f)))
print(f'fontTools vs formula on random cubics: same n {same}, different {diff}, single quadratic {ones}')
for ex in examples:
    print('  fontTools n=%d formula n=%s  |d3|/tol=%.6f  threshold(n_ft)=%.6f threshold(n_f)=%.6f' % ex)
