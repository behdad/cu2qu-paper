#!/usr/bin/env python3
"""The plain-Python reference (cu2qu_geo.py) against the numpy twin (closed.py), and its
equal-parameter n = 1 error (closed form through one cubic equation) against fontTools'
recursive test."""
import sys
sys.path.insert(0, '/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad/ft-old/Lib')
import numpy as np
from common import *
import closed
import cu2qu_geo as G
from fontTools.cu2qu.cu2qu import curve_to_quadratic

rng = np.random.default_rng(5)
sets = {'corpus': load_cubics()}
for name in ['random', 'fontlike', 'cusp', 'loop', 'inflect', 'coincident', 'straight', 'parallel']:
    sets[name] = load_cubics(f'{W}/sets/{name}.bin')

worst = 0.0; nvalid_mismatch = 0; total = 0
for name, c in sets.items():
    idx = rng.choice(len(c), 3000 if name == 'corpus' else 600, replace=False)
    cc = c[idx]
    for n in range(1, 7):
        if n == 1:
            b, v = closed.n1_closed(cc)
        else:
            b, v = closed.ramp_closed(cc, n)
        for j in range(len(cc)):
            pts = [complex(*p) for p in cc[j]]
            s = G.n1_geo_bound(*pts) if n == 1 else G.ramp_geo_bound(*pts, n)
            total += 1
            if (s is None) != (not v[j]):
                # the two may differ only in degenerate corners (d3 = 0, or a leg of length 0)
                nvalid_mismatch += 1
                if nvalid_mismatch <= 5:
                    print('  validity differs:', name, n, cc[j].tolist(), 'scalar', s, 'numpy', b[j], v[j])
                continue
            if s is not None:
                rel = abs(s - b[j]) / max(abs(s), abs(b[j]), 1e-12)
                if rel > worst:
                    worst = rel
                    if rel > 1e-6:
                        print('  differs:', name, n, cc[j].tolist(), 'scalar', s, 'numpy', b[j])
print(f'compared {total} (cubic, n) pairs: validity mismatches {nvalid_mismatch}, largest relative difference of the bounds {worst:.2e}')

# n = 1 equal-parameter error in closed form vs the instrument and vs fontTools' decision
c = sets['corpus']
g = load_gt(W + '/gt_ramp.bin', 16)
bad = 0; worst = 0.0; cnt = 0; dis = 0
for i in range(0, len(c), 7):
    pts = [complex(*p) for p in c[i]]
    q = G.tangent_intersection(*pts)
    if q is None:
        continue
    e = G.n1_ep_error(*pts, q)
    cnt += 1
    if np.isfinite(g[i, 0, 0]):
        worst = max(worst, abs(e - g[i, 0, 0]) / max(e, 1e-9))
    for tol in (1.0, 0.5, 0.1):
        ft1 = len(curve_to_quadratic([tuple(p) for p in c[i]], tol)) == 3
        if ft1 != (e <= tol) and abs(e - tol) > 1e-9:
            dis += 1
print(f'n = 1 equal-parameter error by the cubic formula: {cnt} cubics, largest relative difference from the instrument {worst:.2e}; '
      f'decisions that differ from fontTools (away from ties): {dis}')

# the rule for n, scalar, on a sample of the corpus
H = np.fmax(g[:, :, 1], g[:, :, 2])
for tol in (1.0, 0.5):
    na = first_fit(g[:, :, 0], tol); nc = first_fit(H, tol * (1 + 1e-9))
    ns = []; sel = range(0, len(c), 11)
    for i in sel:
        ns.append(G.quadratics_needed([complex(*p) for p in c[i]], tol))
    ns = np.array(ns); sel = np.array(list(sel))
    print(f'tol {tol}: sample of {len(sel)}: mean n  equal-parameter {na[sel].mean():.4f}  scalar rule {ns.mean():.4f}  truth {nc[sel].mean():.4f};  rule < truth: {(ns < nc[sel]).sum()}')
