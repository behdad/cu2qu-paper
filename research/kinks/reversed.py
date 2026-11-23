"""Single quadratics whose off-curve point is behind an end point: the end
tangents of the cubic meet on the far side of the point, so the quadratic
leaves it in the direction opposite to the cubic's handle.

usage: reversed.py DATA.npz ...
"""

import json
import sys

import numpy as np

total = 0
found = 0
for path in sys.argv[1:]:
    data = np.load(path)
    names = json.load(open(path + ".names.json"))
    pts = data["pts"]
    idx = np.nonzero((data["kind"] == 1) & (data["n"] == 1))[0]
    total += len(idx)
    if not len(idx):
        continue
    p = pts[idx]
    h1 = p[:, :, 1] - p[:, :, 0]
    h2 = p[:, :, 2] - p[:, :, 3]
    g1 = data["qfirst"][idx] - p[:, :, 0]
    g2 = data["qlast"][idx] - p[:, :, 3]
    d1 = (g1 * np.conj(h1)).real
    d2 = (g2 * np.conj(h2)).real
    bad = np.any((d1 < 0) | (d2 < 0), axis=1)
    found += int(np.sum(bad))
    for k in np.nonzero(bad)[0][:5]:
        s = idx[k]
        print("%s: glyph %s segment %d" % (path.split("/")[-1], names[data["gid"][s]], s))
        for m in range(pts.shape[1]):
            print(
                "   master %d: cubic %s  quadratic off-curve %s"
                % (
                    m,
                    " ".join("(%g,%g)" % (z.real, z.imag) for z in pts[s, m]),
                    "(%g,%g)" % (data["qfirst"][s, m].real, data["qfirst"][s, m].imag),
                )
            )
print("single quadratics: %d; with the off-curve point behind an end point in some master: %d" % (total, found))
