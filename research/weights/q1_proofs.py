#!/usr/bin/env python3
"""Q1: exact checks behind the proofs.

A. Lower bound 4/27 whenever ONE end tangent is kept (s_0 = +1), any n >= 1.
   Proof by a chain of point evaluations (see the report).  Here: (i) the two
   inequalities the chain uses, as exact polynomial identities in rationals;
   (ii) a random search that never gets below 4/27.

B. Even n: a dual certificate.  With u_k = 2/3 for even k and 1/3 for odd k and
   sign_k = (-1)^k, there are multipliers lambda_k > 0 such that

       sum_k lambda_k sign_k phi_k(u_k; s)  =  (4/27) sum_k lambda_k

   identically in the free weights s_1..s_{n-2} (s_0 = 1, s_{n-1} = -1).
   Exact rational arithmetic.

C. Odd n: how fast the optimum approaches 4/27; the exact n = 5 optimum.
"""
import sys
from fractions import Fraction as F
import numpy as np
from common import phi_max, alternating


# ---- exact phi in rationals
def phi_exact(s, k, u):
    n = len(s)
    g = lambda j: F(0) if j == 0 or j == n else (s[j - 1] + s[j]) / 8
    h = (1 - s[k]) / 2
    return (1 - u) ** 2 * g(k) + u * u * g(k + 1) + u * (1 - u) * (u - h)


def certificate(n):
    """lambda_0..lambda_{n-1} (Fractions), symmetric, or None."""
    # recurrences from d/ds_j = 0, j = 1..n-2:
    #   j odd : 13 l_j = 4 l_{j-1} + l_{j+1}
    #   j even: 13 l_j = l_{j-1} + 4 l_{j+1}
    def run(l0, l1):
        l = [F(l0), F(l1)]
        for j in range(1, n - 1):
            if j % 2 == 1:
                l.append(13 * l[j] - 4 * l[j - 1])
            else:
                l.append((13 * l[j] - l[j - 1]) / 4)
        return l
    if n == 2:
        return [F(1), F(1)]
    u = run(1, 0); v = run(0, 1)
    # symmetric solution: l_{n-1} = l_0  with l = u + b v
    b = (1 - u[n - 1]) / v[n - 1]
    return [u[k] + b * v[k] for k in range(n)]


def check_certificate(n):
    lam = certificate(n)
    assert all(x > 0 for x in lam), (n, lam)
    assert all(lam[k] == lam[n - 1 - k] for k in range(n))
    uk = [F(2, 3) if k % 2 == 0 else F(1, 3) for k in range(n)]
    sg = [1 if k % 2 == 0 else -1 for k in range(n)]
    def combo(s):
        return sum(lam[k] * sg[k] * phi_exact(s, k, uk[k]) for k in range(n))
    base = [F(1 if k % 2 == 0 else -1) for k in range(n)]
    c0 = combo(base)
    assert c0 == F(4, 27) * sum(lam), (n, c0)
    # the combination is affine in s, so independence of each free s_j is
    # checked by moving one weight at a time
    for j in range(1, n - 1):
        s = list(base); s[j] += F(7, 3)
        assert combo(s) == c0, (n, j)
    return lam


if __name__ == '__main__':
    print('A. chain inequalities, exact:')
    # even piece: phi_k(2/3) - 4/27 - 4 g_{k+1}/9 = g_k/9 + (s_k - 1)/9   (>= 0 when g_k >= 0, s_k >= 1)
    # odd piece : phi_k(1/3) + 4/27 - g_{k+1}/9  = 4 g_k/9 + (s_k + 1)/9  (<= 0 when g_k <= 0, s_k <= -1)
    import random
    random.seed(5)
    for _ in range(2000):
        n = random.randint(3, 9)
        s = [F(random.randint(-40, 40), random.randint(1, 9)) for _ in range(n)]
        k = random.randint(1, n - 2)
        gk = (s[k - 1] + s[k]) / 8; gk1 = (s[k] + s[k + 1]) / 8
        assert phi_exact(s, k, F(2, 3)) - F(4, 27) - 4 * gk1 / 9 == gk / 9 + (s[k] - 1) / 9
        assert phi_exact(s, k, F(1, 3)) + F(4, 27) - gk1 / 9 == 4 * gk / 9 + (s[k] + 1) / 9
    print('   phi_k(2/3) = 4/27 + 4 g_{k+1}/9 + g_k/9 + (s_k - 1)/9      identity holds')
    print('   phi_k(1/3) = -4/27 + g_{k+1}/9 + 4 g_k/9 + (s_k + 1)/9     identity holds')
    rng = np.random.default_rng(7)
    best = 1e9
    for _ in range(200000):
        n = int(rng.integers(1, 9))
        s = rng.normal(0, 1.5, size=n) * rng.choice([0.2, 1, 3])
        # start near alternating half the time
        if rng.random() < 0.5:
            s = alternating(n) * (1 + np.abs(rng.normal(0, 0.3, size=n)))
        s[0] = 1.0
        best = min(best, phi_max(s))
    print(f'   random search, s_0 = 1, n = 1..8, 200000 trials: least max|phi| found {best:.9f}  (4/27 = {4/27:.9f})')

    print()
    print('B. dual certificate for even n (exact rationals):')
    for n in (2, 4, 6, 8, 10, 12, 20, 40, 100):
        lam = check_certificate(n)
        tot = sum(lam)
        txt = ', '.join(str(x / lam[n // 2]) for x in lam[:min(n, 6)])
        print(f'   n = {n:3d}: all lambda > 0, identity exact;  lambda/lambda_mid = {txt}{" ..." if n > 6 else ""}')
    print('   transfer matrix over a pair [[-4, 13], [-13, 42]], eigenvalues 19 +- 6 sqrt(10) =',
          19 + 6 * 10 ** 0.5, 19 - 6 * 10 ** 0.5)

    print()
    print('C. odd n: optimum minus 4/27, and the ratio between n and n + 4')
    from q1_lp import solve
    vals = {}
    for n in range(3, 24, 2):
        t, s = solve(n, 'both', antisym=True)
        vals[n] = t - 4 / 27
    for n in range(3, 24, 2):
        r = f'{vals[n - 4] / vals[n]:.3f}' if n - 4 in vals and vals[n] > 1e-11 else '   -'
        print(f'   n = {n:2d}: Phi - 4/27 = {vals[n]:.3e}   ratio to n-4: {r}')
    print('   (19 + 6 sqrt(10) = %.4f)' % (19 + 6 * 10 ** 0.5))
