#!/usr/bin/env python3
"""The rule for n in the plain-Python reference (cu2qu_geo.quadratics_needed) against the
numpy pipeline that produced the tables, cubic by cubic, on every 5th corpus cubic."""
import numpy as np
from common import *
from closed import *
import cu2qu_geo as G

c = load_cubics()
NMAX = 16
g = load_gt(W + '/gt_ramp.bin', NMAX)
ep = g[:, :, 0]
sel = np.arange(0, len(c), 5)
cc = c[sel]
N = len(cc)
d3 = np.hypot(*d3_of(cc).T)
B = np.full((N, NMAX), np.inf)
b, v = n1_closed(cc)
B[:, 0] = np.where(v, np.minimum(b, ep[sel, 0]), ep[sel, 0])
for n in range(2, NMAX + 1):
    b, v = ramp_closed(cc, n)
    e = d3 * Phi_ramp(n) / n ** 3
    B[:, n - 1] = np.where(v, np.minimum(b, e), e)
for tol in (1.0, 0.5, 0.1):
    nb = first_fit(B, tol)
    ns = np.array([G.quadratics_needed([complex(*p) for p in x], tol, max_n=NMAX) or 0 for x in cc])
    diff = np.where(nb != ns)[0]
    print(f'tol {tol}: {N} cubics, mean n numpy {nb.mean():.4f} scalar {ns.mean():.4f}, cubics where they differ: {len(diff)}')
    for j in diff[:5]:
        print('    ', cc[j].tolist(), 'numpy', nb[j], 'scalar', ns[j], 'B', B[j, :4])
