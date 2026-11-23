#!/usr/bin/env python3
"""Average quadratics per cubic for the numerical yardsticks against the
closed-form schemes, on the cubics the yardstick was run on (the whole corpus,
or a sample: pass the sample tag, e.g. `_sample2000`).

Also: how far each closed-form scheme's error is from the yardstick at equal n.
"""
import json, os, sys
import numpy as np
from common import CORPUS
from corpus_stats import least_n, TOLS, BIG

tag = sys.argv[1] if len(sys.argv) > 1 else ''
T = np.load(os.path.join(BIG, f'corpus_yard_tan{tag}.npy'))
Fr = np.load(os.path.join(BIG, f'corpus_yard_free{tag}.npy'))
idx = np.load(os.path.join(BIG, f'corpus_yard_idx{tag}.npy')) if tag else np.arange(len(T))
N = len(idx)
ft_n = np.load(os.path.join(BIG, 'corpus_ft_n.npy'))[idx]
names = [('ramp', 'ramp'), ('optu', 'optu'), ('Tuni', 'T equal-parameter knots'), ('Tgam', 'T gamma knots'),
         ('Taff', 'T affine knots'), ('relax50', 'relax50'), ('zero', 'zero')]
H = {k: np.load(os.path.join(BIG, f'corpus_{k}.npy'))[idx] for k, _ in names}

def filled(A):
    """n not computed because an earlier n was already below the stop: treat as -inf (fits)."""
    A = A.copy()
    done = np.zeros(len(A), dtype=bool)
    for j in range(A.shape[1]):
        col = A[:, j]
        col[np.isnan(col) & done] = 0.0
        done |= (col <= 0.1)
        A[:, j] = col
    A[np.isnan(A)] = np.inf
    return A

Tf = filled(T); Ff = filled(Fr)
print(f'{N} cubics.  Average quadratics per cubic (n chosen by measured Hausdorff distance):')
print(f'  {"scheme":40s} ' + ''.join(f'tol {t:<8}' for t in TOLS))
print(f'  {"cu2qu as implemented":40s} ' + ''.join(f'{ft_n[:, j].mean():<12.4f}' for j in range(len(TOLS))))
for k, label in names:
    print(f'  {label:40s} ' + ''.join(f'{least_n(H[k], t).mean():<12.4f}' for t in TOLS))
for label, A in (('TAN*  (numerical optimum, tangents kept)', Tf), ('FREE* (numerical optimum, tangents free)', Ff)):
    print(f'  {label:40s} ' + ''.join(f'{least_n(A, t).mean():<12.4f}' for t in TOLS))
print()
print('Distribution of n, percent (n = 1..6, 7+):')
for label, A in (('TAN*', Tf), ('FREE*', Ff)):
    for t in TOLS:
        cnt = np.bincount(np.minimum(least_n(A, t), 7), minlength=8)
        print(f'  {label:6s} tol {t:<4}: ' + ' '.join(f'{100 * cnt[j] / N:5.1f}' for j in range(1, 8)))
print()
print('Error of a scheme divided by the yardstick error at the same n (cubics where both were computed')
print('and the yardstick error is above 0.02 units): median, and 90th percentile')
for k, label in names:
    yard = Tf if k in ('ramp', 'optu', 'Tuni', 'Tgam', 'Taff') else Ff
    line = []
    for n in (1, 2, 3, 4):
        a = H[k][:, n - 1]; b = yard[:, n - 1]
        ok = np.isfinite(a) & np.isfinite(b) & (b > 0.02)
        if ok.sum() < 20:
            line.append('      -      ')
        else:
            r = a[ok] / b[ok]
            line.append(f'{np.median(r):5.2f} {np.quantile(r, 0.9):6.2f} ')
    print(f'  {label:26s} vs {"TAN* " if yard is Tf else "FREE*"}  n=1: {line[0]}  n=2: {line[1]}  n=3: {line[2]}  n=4: {line[3]}')
