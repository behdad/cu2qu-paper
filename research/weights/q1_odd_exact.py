#!/usr/bin/env python3
"""Q1: the exact optimum for n = 5, and 30-digit values for odd n up to 13.

Odd n = 2m + 1, weights antisymmetric (s_k = -s_{n-1-k}, s_m = 0).  At the
optimum pieces 0..m-1 each touch +-t once, with alternating signs, and the
middle piece is slack.  For n = 5:  s = (1, a, 0, -a, -1),

    t = (4/27) ((9 + a) / 8)^3              (piece 0)
    min_u phi_1(u; a) = -t                  (piece 1)

so a is a root of  disc_u (phi_1(u; a) + t(a)) = 0.

Run with the venv python (needs sympy, mpmath):
  /home/behdad/hb/cu2qu-work/research-weights/venv/bin/python q1_odd_exact.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import mpmath as mp

a, u = sp.symbols('a u', real=True)
t = sp.Rational(4, 27) * ((9 + a) / 8) ** 3
g1 = (1 + a) / 8
g2 = a / 8
h1 = (1 - a) / 2
phi1 = (1 - u) ** 2 * g1 + u ** 2 * g2 + u * (1 - u) * (u - h1)
D = sp.discriminant(sp.expand(phi1 + t), u)
num = sp.numer(sp.together(D))
facs = sp.factor_list(sp.expand(num))
print('disc_u(phi_1 + t), numerator factors:')
roots = []
for f, mult in facs[1]:
    print('   (', sp.Poly(f, a).as_expr(), ') ^', mult)
    for r in sp.Poly(f, a).nroots(n=30):
        if abs(sp.im(r)) < 1e-25:
            roots.append((sp.re(r), f))
for r, f in sorted(roots, key=lambda x: x[0]):
    tv = sp.N(t.subs(a, r), 30)
    print(f'   real root a = {sp.N(r, 25)}   t = {tv}')
best = min(roots, key=lambda x: abs(x[0] + sp.Rational(9334, 10000)))
r, f = best
print()
print('n = 5 optimum: a = s_1 =', sp.N(r, 30))
print('               Phi = t =', sp.N(t.subs(a, r), 30))
print('               a is a root of', sp.Poly(f, a).as_expr(), '= 0')
T = sp.symbols('T')
res = sp.resultant(sp.Poly(f, a).as_expr(), sp.numer(sp.together(T - t)), a)
for ff, mult in sp.factor_list(res)[1]:
    for rr in sp.Poly(ff, T).nroots(n=30):
        if abs(sp.im(rr)) < 1e-25 and abs(sp.re(rr) - t.subs(a, r)) < 1e-20:
            print('               Phi is a root of', sp.Poly(ff, T).as_expr(), '= 0')

# ---- high precision for n = 5..13: unknowns s_1..s_{m-1}, t
mp.mp.dps = 40

def piece_ext(gk, gk1, hk, want_max):
    """extreme value of (1-u)^2 gk + u^2 gk1 + u (1-u) (u - hk) at its critical points"""
    c0, c1, c2 = gk, -(2 * gk + hk), 1 + gk + gk1 + hk
    disc = (2 * c2) ** 2 + 12 * c1
    rt = mp.sqrt(disc)
    us = [(2 * c2 - rt) / 6, (2 * c2 + rt) / 6]
    vals = [c0 + x * (c1 + x * (c2 - x)) for x in us]
    return max(vals) if want_max else min(vals)

def residual(n, x):
    m = (n - 1) // 2
    s = [mp.mpf(1)] + list(x[:m - 1]) + [mp.mpf(0)]
    tt = x[m - 1]
    out = []
    for k in range(m):
        gk = 0 if k == 0 else (s[k - 1] + s[k]) / 8
        gk1 = (s[k] + s[k + 1]) / 8
        hk = (1 - s[k]) / 2
        if k % 2 == 0:
            out.append(piece_ext(gk, gk1, hk, True) - tt)
        else:
            out.append(piece_ext(gk, gk1, hk, False) + tt)
    return out

from q1_lp import solve
print()
print('odd n, 30 digits (Newton on the equioscillation system, started from the LP):')
for n in (5, 7, 9, 11, 13):
    m = (n - 1) // 2
    t0, s0 = solve(n, 'both', antisym=True)
    x0 = [mp.mpf(float(v)) for v in s0[1:m]] + [mp.mpf(float(t0))]
    sol = mp.findroot(lambda *x: residual(n, x), x0, tol=1e-32, maxsteps=50)
    xs = list(sol)
    print(f'   n = {n:2d}: Phi = {mp.nstr(xs[-1], 30)}   n^3 / Phi = {mp.nstr(n ** 3 / xs[-1], 12)}')
    print('           s_1..s_{m-1} = ' + ', '.join(mp.nstr(v, 20) for v in xs[:-1]))
