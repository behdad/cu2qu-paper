"""End-to-end check of stage 2 against a real build.

usage: endtoend.py DATA.npz GLYPH [GLYPH ...]

The cubic masters of the given glyphs (two masters, one axis) are written to
two UFOs in memory, built into a TrueType variable font with ufo2ft (which
converts with fontTools.cu2qu and rounds to integers), and the font's outlines
are read back at locations along the axis.  The angle between the two handles
at every joint is then measured in the font and compared with what stage 2
computes for its rounded variant at the same locations.
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stage2

DEG = 180.0 / math.pi


def contours_of(data, gi):
    """Per master, the contours of glyph gi as lists of segments."""
    segs = np.nonzero(data["gid"] == gi)[0]
    prev = data["prev"]
    out = []
    cur = []
    for s in segs:
        if cur and prev[s] != s - 1:
            out.append(cur)
            cur = []
        cur.append(s)
    if cur:
        out.append(cur)
    return out


def main():
    import ufoLib2
    from fontTools.designspaceLib import (
        AxisDescriptor,
        DesignSpaceDocument,
        SourceDescriptor,
    )
    from fontTools.pens.recordingPen import RecordingPen
    import ufo2ft

    path = sys.argv[1]
    glyphs = sys.argv[2:]
    data, meta, names = stage2.load(path)
    pts = data["pts"]
    kind = data["kind"]
    M = pts.shape[1]
    assert M == 2 and len(meta["locations"]) == 2
    res = stage2.analyze(data, meta)
    joint_of = {int(j): k for k, j in enumerate(res["j"])}

    fonts = []
    for m in range(M):
        font = ufoLib2.Font()
        font.info.unitsPerEm = meta["upem"]
        font.info.ascender = 880
        font.info.descender = -120
        font.info.xHeight = 500
        font.info.capHeight = 700
        font.info.familyName = "KinkTest"
        font.info.styleName = "Master%d" % m
        g = font.newGlyph(".notdef")
        g.width = 500
        for name in glyphs:
            gi = names.index(name)
            g = font.newGlyph(name)
            g.width = 1000
            pen = g.getPen()
            for contour in contours_of(data, gi):
                p0 = pts[contour[0], m, 0]
                pen.moveTo((p0.real, p0.imag))
                for s in contour:
                    if kind[s] == 0:
                        p = pts[s, m, 1]
                        pen.lineTo((p.real, p.imag))
                    else:
                        pen.curveTo(*[(z.real, z.imag) for z in pts[s, m, 1:]])
                pen.closePath()
        fonts.append(font)

    doc = DesignSpaceDocument()
    axis = AxisDescriptor()
    axis.name = "Weight"
    axis.tag = "wght"
    axis.minimum, axis.default, axis.maximum = 0, 0, 1000
    doc.addAxis(axis)
    for m, font in enumerate(fonts):
        src = SourceDescriptor()
        src.font = font
        src.location = {"Weight": 1000 * m}
        src.familyName = "KinkTest"
        src.styleName = "Master%d" % m
        doc.addSource(src)
    vf = ufo2ft.compileVariableTTF(doc, reverseDirection=False)

    worst = 0.0
    count = 0
    biggest = (0.0, None)
    for li, loc in enumerate(res["samples"]):
        t = list(loc.values())[0]
        gs = vf.getGlyphSet(location={"wght": 1000 * t})
        w = res["W"][li]
        for name in glyphs:
            gi = names.index(name)
            pen = RecordingPen()
            gs[name].draw(pen)
            # The font's contours as point lists with on/off marks.
            fcontours = []
            for op, args in pen.value:
                if op == "moveTo":
                    cur = [(complex(*args[0]), True)]
                elif op == "lineTo":
                    cur.append((complex(*args[0]), True))
                elif op == "qCurveTo":
                    for a in args[:-1]:
                        cur.append((complex(*a), False))
                    cur.append((complex(*args[-1]), True))
                elif op == "closePath":
                    if cur[-1][0] == cur[0][0] and cur[-1][1]:
                        cur.pop()
                    fcontours.append(cur)
            for contour, fc in zip(contours_of(data, gi), fcontours):
                for s in contour:
                    k = joint_of[int(s)]
                    # The joint in the font: masters rounded, then interpolated.
                    on = sum(w[m] * stage2.oround(pts[s, m, 0]) for m in range(M))
                    idx = [i for i, (z, o) in enumerate(fc) if o and abs(z - on) < 1e-6]
                    if len(idx) != 1:
                        continue
                    i = idx[0]
                    before = fc[i - 1][0]
                    after = fc[(i + 1) % len(fc)][0]
                    inh = on - before
                    outh = after - on
                    if inh == 0 or outh == 0:
                        continue
                    font_kink = math.atan2(
                        (outh * inh.conjugate()).imag, (outh * inh.conjugate()).real
                    ) * DEG
                    # Stage 2's rounded variant at this location.
                    mine = one_kink(data, int(s), w)
                    if mine is None:
                        continue
                    d = abs(stage2.wrap(np.array([font_kink - mine]))[0])
                    worst = max(worst, d)
                    count += 1
                    if abs(font_kink) > biggest[0]:
                        biggest = (abs(font_kink), (name, int(s), t, font_kink, mine))
    print("%d joint samples compared; largest difference between the font and stage 2: %.3g degrees" % (count, worst))
    print("largest kink seen in the font: %.4f degrees at %s" % (biggest[0], biggest[1]))


def one_kink(data, s, w):
    """The kink of the rounded quadratic variant at the start of segment s,
    for master weights w, computed from the stored first and last off-curve
    points."""
    pts = data["pts"]
    kind = data["kind"]
    p = int(data["prev"][s])
    M = pts.shape[1]
    inh = outh = 0j
    for m in range(M):
        on = stage2.oround(pts[s, m, 0])
        if kind[s] == 1:
            o = stage2.oround(data["qfirst"][s, m]) - on
        else:
            o = stage2.oround(pts[s, m, 1]) - on
        if kind[p] == 1:
            i = on - stage2.oround(data["qlast"][p, m])
        else:
            i = on - stage2.oround(pts[p, m, 0])
        inh += w[m] * i
        outh += w[m] * o
    if inh == 0 or outh == 0:
        return None
    z = outh * inh.conjugate()
    return math.atan2(z.imag, z.real) * DEG


if __name__ == "__main__":
    main()
