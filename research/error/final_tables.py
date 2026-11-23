#!/usr/bin/env python3
"""The tables of the report, corpus part.

  (a) equal-parameter error (what cu2qu tests today)
  (b) the closed-form geometric bound, with (a) as the fallback where it is not valid
  (c) the Hausdorff distance measured by gt.c

Usage: final_tables.py [table...]   (1 = ratios per n, 2 = quadratics per cubic, 3 = n histogram)
"""
import sys
import numpy as np
from common import *
from closed import *
from geo import pointwise_dense, weights
import n1_identity

which = set(sys.argv[1:]) or {'1', '2', '3'}
c = load_cubics()
NMAX = 16
g = load_gt(W + '/gt_ramp.bin', NMAX)
H = np.fmax(g[:, :, 1], g[:, :, 2])
ep = g[:, :, 0]
N = len(c)
d3 = np.hypot(*d3_of(c).T)
EPS = 1e-9

def build(variant):
    """(N, NMAX) array of the bound used by a variant."""
    B = np.full((N, NMAX), np.inf); V = np.zeros((N, NMAX), bool)
    for n in range(1, NMAX + 1):
        epn = ep[:, 0] if n == 1 else d3 * Phi_ramp(n) / n ** 3
        if n == 1:
            if variant == 'geo_n>=2 only':
                b, v = np.full(N, np.inf), np.zeros(N, bool)
            else:
                b, v = n1_closed(c, refined=(variant in ('final', 'geo_n=1 only')))
        else:
            if variant == 'geo_n=1 only':
                b, v = np.full(N, np.inf), np.zeros(N, bool)
            elif variant == 'product':
                b, v = ramp_closed(c, n, engine='prod')
            elif variant == 'simple':
                b, v = ramp_closed(c, n, refine=0)
            else:
                b, v = ramp_closed(c, n, refine=1)
        B[:, n - 1] = np.where(v, np.minimum(b, epn), epn); V[:, n - 1] = v
    return B, V

def build_dense():
    """The pointwise bounds maximised by dense sampling instead of in closed form (not a formula:
    it shows what an ideal maximisation of the same pointwise bound would give)."""
    B = np.full((N, NMAX), np.inf)
    sd, first, bad, (l1, l2, p, r) = n1_identity.sampson_dense(c, T=1001)
    v = np.isfinite(l1) & (l1 > 0) & (l1 < 2) & (l2 > 0) & (l2 < 2) & ~bad & np.isfinite(sd)
    B[:, 0] = np.where(v, np.minimum(sd, ep[:, 0]), ep[:, 0])
    for n in range(2, NMAX + 1):
        epn = d3 * Phi_ramp(n) / n ** 3
        m = np.isfinite(H[:, n - 1])
        b = np.full(N, np.inf)
        if m.any():
            f, e, bad = pointwise_dense(c[m], weights('ramp', n), U=257)
            b[m] = np.where(bad, np.inf, e)
        B[:, n - 1] = np.minimum(b, epn)
    return B

Bf, Vf = build('final')

if '1' in which:
    print('Table 1.  Corpus, ramp.  Ratio of each bound to the measured Hausdorff distance, per n.')
    print('          (cubics for which n was computed and the distance is > 1e-6)')
    print(' n   cubics  (b) valid   (b) under |  (a)/(c):   5%    25%    50%    75%    95%  |  (b)/(c):  min       5%     25%     50%     75%     95%     99%')
    for n in range(1, 11):
        m = np.isfinite(H[:, n - 1]) & (H[:, n - 1] > 1e-6) & np.isfinite(ep[:, n - 1])
        if not m.any():
            continue
        ra = ep[m, n - 1] / H[m, n - 1]; rb = Bf[m, n - 1] / H[m, n - 1]
        under = (Bf[m, n - 1] < H[m, n - 1] * (1 - EPS)).sum()
        qa = np.percentile(ra, [5, 25, 50, 75, 95]); qb = np.percentile(rb, [5, 25, 50, 75, 95, 99])
        print(f'{n:2d}  {m.sum():6d}   {100 * Vf[m, n - 1].mean():5.1f}%    {under:5d}     |          ' + ' '.join(f'{x:6.2f}' for x in qa) +
              f'  |        {rb.min():.6f} ' + ' '.join(f'{x:7.4f}' for x in qb))
    # all n together, at the n each cubic actually needs
    for tol in (1.0, 0.5, 0.1):
        nc = first_fit(H, tol * (1 + EPS))
        ok = nc > 0
        i = np.arange(N)[ok]; k = nc[ok] - 1
        h = H[i, k]; m = h > 1e-6
        ra = ep[i, k][m] / h[m]; rb = Bf[i, k][m] / h[m]
        print(f'   at the true least n for tol {tol}: (a)/(c) median {np.median(ra):.3f}, 90% {np.percentile(ra, 90):.3f};  (b)/(c) median {np.median(rb):.4f}, 90% {np.percentile(rb, 90):.4f}, 99% {np.percentile(rb, 99):.3f}')
    print()

if '2' in which:
    variants = [('product', 'n=1: M*F;  n>=2: Phi_k * max sin for every piece'),
                ('simple', 'n=1: M*F;  n>=2: first/last piece maximised, one tangent point'),
                ('final', 'n=1 and first/last piece with linear majorants, tangent point refined once'),
                ('geo_n=1 only', 'geometric test for n = 1 only (n >= 2 as today)'),
                ('geo_n>=2 only', 'geometric test for n >= 2 only (n = 1 as today)')]
    Bs = {'final': Bf}
    for v, _ in variants:
        if v not in Bs:
            Bs[v] = build(v)[0]
    Bd = build_dense()
    print('Table 2.  Corpus (77,533 cubics), ramp.  Mean quadratics per cubic; in brackets the saving against (a).')
    print(f'{"measure":78s}' + ''.join(f'   tol {t:<4}        ' for t in (1.0, 0.5, 0.1)))
    rows = [('(a) equal-parameter (cu2qu today)', ep)]
    rows += [(f'(b) closed form, {v}: {txt}'[:78], Bs[v]) for v, txt in variants]
    rows += [('(b*) same pointwise bounds maximised by dense sampling (not a formula)', Bd), ('(c) true least n (measured Hausdorff distance)', H)]
    base = {}; notes = []
    for name, arr in rows:
        line = f'{name:78s}'
        for tol in (1.0, 0.5, 0.1):
            na = first_fit(ep, tol)
            ok = na > 0
            nn = first_fit(arr, tol * (1 + EPS) if arr is H else tol)
            mean = nn[ok].mean()
            base.setdefault(tol, mean)
            line += f'   {mean:.4f} ({100 * (1 - mean / base[tol]):4.1f}%)'
            if arr is not H and arr is not ep:
                nc = first_fit(H, tol * (1 + EPS))
                fewer = (nn[ok] < nc[ok]).sum()
                if arr is Bd:
                    notes.append(f'(b*) picks fewer than (c) for {fewer} cubics at tol {tol} (sampling of the maximum)')
                else:
                    assert fewer == 0, (name, tol, fewer)
        print(line)
    print('   (every (b) row was checked to never pick fewer quadratics than (c); ' + '; '.join(notes) + ')')
    print()

if '3' in which:
    print('Table 3.  Corpus, ramp.  How many cubics take n quadratics.')
    for tol in (1.0, 0.5, 0.1):
        na = first_fit(ep, tol); nb = first_fit(Bf, tol); nc = first_fit(H, tol * (1 + EPS))
        ok = na > 0
        print(f'  tol {tol}:     n =' + ''.join(f'{k:7d}' for k in range(1, 11)) + '    more')
        for name, nn in (('(a)', na), ('(b)', nb), ('(c)', nc)):
            print(f'       {name}      ' + ''.join(f'{(nn[ok] == k).sum():7d}' for k in range(1, 11)) + f'{(nn[ok] > 10).sum():7d}')
        print(f'       cubics where (b) needs more than (c): {(nb[ok] > nc[ok]).sum()} ({100 * (nb[ok] > nc[ok]).mean():.2f}%);  fewer: {(nb[ok] < nc[ok]).sum()}')
    print()
