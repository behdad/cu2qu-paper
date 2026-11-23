"""Shared code for the cu2qu kink study.

Stage 1 (this file, `extract_*` and `convert`): get point-compatible cubic
masters out of a CFF2 variable font, a designspace of UFOs or a .glyphs file,
as flat numpy arrays, and convert every cubic across masters with
`fontTools.cu2qu.cu2qu.curves_to_quadratic`.

The flat layout, for S segments and M masters:

  pts[S, M, 4]  complex128  control points; a line has p0, p1, nan, nan
  kind[S]       uint8       0 line, 1 cubic
  prev[S]       int64       index of the segment that ends where this starts
                            (-1 at the start of an open contour)
  flag[S]       int8        1 if the point where the segment starts is marked
                            smooth in every master, 0 if not, -1 if the source
                            has no such marks
  gid[S]        int32       index into names
  n[S]          int16       quadratics per cubic that curves_to_quadratic gave
                            (0 for a line)
  qfirst[S, M], qlast[S, M] complex128  first and last off-curve point of the
                            spline (nan for a line)
"""

import math
import sys

import numpy as np

NAN = float("nan")


# ---------------------------------------------------------------------------
# extraction
# ---------------------------------------------------------------------------


class _Flat:
    def __init__(self, M):
        self.M = M
        self.pts = []
        self.kind = []
        self.prev = []
        self.flag = []
        self.gid = []
        self.names = []
        self.skipped = {}

    def skip(self, why):
        self.skipped[why] = self.skipped.get(why, 0) + 1

    def add_glyph(self, name, contours):
        """contours: list of (closed, segs, flags); segs[i][m] is a tuple of 2
        or 4 complex; flags[i] is the mark for the point where segs[i] starts.
        """
        g = len(self.names)
        self.names.append(name)
        for closed, segs, flags in contours:
            if not segs:
                continue
            base = len(self.kind)
            count = len(segs)
            for i, seg in enumerate(segs):
                k = len(seg[0])
                row = []
                for m in range(self.M):
                    p = seg[m]
                    if k == 2:
                        row.append((p[0], p[1], complex(NAN, NAN), complex(NAN, NAN)))
                    else:
                        row.append(tuple(p))
                self.pts.append(row)
                self.kind.append(0 if k == 2 else 1)
                if i == 0:
                    self.prev.append(base + count - 1 if closed else -1)
                else:
                    self.prev.append(base + i - 1)
                self.flag.append(flags[i])
                self.gid.append(g)

    def arrays(self):
        return dict(
            pts=np.array(self.pts, dtype=np.complex128).reshape(-1, self.M, 4),
            kind=np.array(self.kind, dtype=np.uint8),
            prev=np.array(self.prev, dtype=np.int64),
            flag=np.array(self.flag, dtype=np.int8),
            gid=np.array(self.gid, dtype=np.int32),
            names=np.array(self.names, dtype=object),
        )


def _contours_from_recordings(recs):
    """recs: per master, a RecordingPen.value. Returns contours or a string
    saying why the glyph is skipped."""
    M = len(recs)
    ops = [op for op, _ in recs[0]]
    for r in recs[1:]:
        if [op for op, _ in r] != ops:
            return "incompatible operators"
    contours = []
    segs = None
    start = cur = None
    for i, op in enumerate(ops):
        args = [r[i][1] for r in recs]
        if op == "moveTo":
            if segs:
                return "moveTo inside a contour"
            segs = []
            start = [complex(*a[0]) for a in args]
            cur = list(start)
        elif op == "lineTo":
            end = [complex(*a[0]) for a in args]
            segs.append([(cur[m], end[m]) for m in range(M)])
            cur = end
        elif op == "curveTo":
            if any(len(a) != 3 for a in args):
                return "curveTo without three points"
            pts = [[complex(*p) for p in a] for a in args]
            segs.append(
                [(cur[m], pts[m][0], pts[m][1], pts[m][2]) for m in range(M)]
            )
            cur = [p[2] for p in pts]
        elif op in ("closePath", "endPath"):
            closed = op == "closePath"
            if closed and any(cur[m] != start[m] for m in range(M)):
                segs.append([(cur[m], start[m]) for m in range(M)])
            # A line of no length in every master is no segment.
            segs = [
                s
                for s in segs
                if not (len(s[0]) == 2 and all(p[0] == p[1] for p in s))
            ]
            contours.append((closed, segs, [-1] * len(segs)))
            segs = None
        else:
            return "operator " + op
    if segs:
        return "contour not ended"
    return contours


def extract_cff2(path, font_number, glyph_range=None, locations=None):
    """Masters of a CFF2 variable font at the peaks of its regions.

    Returns (arrays, master_locations, upem, skipped)."""
    from fontTools.pens.recordingPen import RecordingPen
    from fontTools.ttLib import TTFont

    font = TTFont(path, fontNumber=font_number) if font_number >= 0 else TTFont(path)
    tags = [a.axisTag for a in font["fvar"].axes]
    if locations is None:
        store = font["CFF2"].cff.topDictIndex[0].VarStore.otVarStore
        peaks = set()
        for region in store.VarRegionList.Region:
            peaks.add(tuple(a.PeakCoord for a in region.VarRegionAxis))
        peaks.add(tuple(0.0 for _ in tags))
        locations = [dict(zip(tags, p)) for p in sorted(peaks)]
    glyphsets = [font.getGlyphSet(location=loc, normalized=True) for loc in locations]
    order = font.getGlyphOrder()
    if glyph_range is not None:
        order = order[glyph_range[0] : glyph_range[1]]
    flat = _Flat(len(locations))
    for name in order:
        recs = []
        for gs in glyphsets:
            pen = RecordingPen()
            gs[name].draw(pen)
            recs.append(pen.value)
        contours = _contours_from_recordings(recs)
        if isinstance(contours, str):
            flat.skip(contours)
            continue
        flat.add_glyph(name, contours)
    return flat.arrays(), locations, font["head"].unitsPerEm, flat.skipped


def _otround(z):
    return complex(math.floor(z.real + 0.5), math.floor(z.imag + 0.5))


def extract_ttvf(path, font_number, glyph_range=None, locations=None, rounded=True):
    """Cubic masters made from a TrueType variable font: every quadratic of
    its outlines at the ends (and the default) of its axes is raised to a
    cubic, and the points are rounded to integers if `rounded`, the way a
    cubic source drawn over TrueType outlines would hold them.

    Returns (arrays, master_locations, upem, skipped)."""
    import itertools

    from fontTools.pens.basePen import decomposeQuadraticSegment
    from fontTools.pens.recordingPen import RecordingPen
    from fontTools.ttLib import TTFont

    font = TTFont(path, fontNumber=font_number) if font_number >= 0 else TTFont(path)
    axes = font["fvar"].axes
    if locations is None:
        per_axis = []
        for a in axes:
            vals = [0.0]
            if a.minValue < a.defaultValue:
                vals.insert(0, -1.0)
            if a.maxValue > a.defaultValue:
                vals.append(1.0)
            per_axis.append(vals)
        if len(axes) > 2:
            # Only the default and the ends of each axis alone.
            combos = [tuple(0.0 for _ in axes)]
            for i, vals in enumerate(per_axis):
                for v in vals:
                    if v:
                        combos.append(tuple(v if k == i else 0.0 for k in range(len(axes))))
        else:
            combos = list(itertools.product(*per_axis))
        locations = [dict(zip([a.axisTag for a in axes], c)) for c in combos]
    glyphsets = [font.getGlyphSet(location=loc, normalized=True) for loc in locations]
    order = font.getGlyphOrder()
    if glyph_range is not None:
        order = order[glyph_range[0] : glyph_range[1]]
    M = len(locations)
    flat = _Flat(M)
    fix = _otround if rounded else (lambda z: z)
    for name in order:
        recs = []
        for gs in glyphsets:
            pen = RecordingPen()
            gs[name].draw(pen)
            recs.append(pen.value)
        ops = [op for op, _ in recs[0]]
        if "addComponent" in ops:
            flat.skip("composite")
            continue
        if any([op for op, _ in r] != ops for r in recs[1:]):
            flat.skip("incompatible operators")
            continue
        # Rewrite as lines and cubics.
        out = [[] for _ in range(M)]
        bad = None
        for i, op in enumerate(ops):
            for m in range(M):
                args = recs[m][i][1]
                if op == "qCurveTo":
                    if args[-1] is None:
                        bad = "no on-curve point"
                        break
                    cur = out[m][-1][1][-1]
                    for q, p2 in decomposeQuadraticSegment(args):
                        p0 = complex(*cur)
                        qq = complex(*q)
                        pp = complex(*p2)
                        c1 = fix(p0 + (qq - p0) * (2.0 / 3.0))
                        c2 = fix(pp + (qq - pp) * (2.0 / 3.0))
                        pp = fix(pp)
                        cur = (pp.real, pp.imag)
                        out[m].append(
                            ("curveTo", ((c1.real, c1.imag), (c2.real, c2.imag), cur))
                        )
                elif op in ("moveTo", "lineTo"):
                    p = fix(complex(*args[0]))
                    out[m].append((op, ((p.real, p.imag),)))
                else:
                    out[m].append((op, args))
            if bad:
                break
        if bad:
            flat.skip(bad)
            continue
        contours = _contours_from_recordings(out)
        if isinstance(contours, str):
            flat.skip(contours)
            continue
        if contours:
            flat.add_glyph(name, contours)
    return flat.arrays(), locations, font["head"].unitsPerEm, flat.skipped


def _ufo_contour(contour):
    """A UFO contour as (closed, [(kind, points, smooth)]) where points are the
    control points after the segment's start and smooth is the mark of the
    segment's end point; or a string saying why not."""
    points = list(contour.points)
    if not points:
        return "empty contour"
    closed = points[0].type != "move"
    if closed:
        for i, p in enumerate(points):
            if p.type is not None:
                break
        else:
            return "no on-curve point"
        points = points[i + 1 :] + points[: i + 1]
        start = points[-1]
    else:
        start = points[0]
        points = points[1:]
    out = []
    offs = []
    for p in points:
        if p.type is None:
            offs.append(complex(p.x, p.y))
            continue
        if p.type == "line" or (p.type == "curve" and not offs):
            if offs:
                return "line with off-curve points"
            out.append((0, (complex(p.x, p.y),), bool(p.smooth)))
        elif p.type == "curve":
            if len(offs) != 2:
                return "curve without two off-curve points"
            out.append((1, (offs[0], offs[1], complex(p.x, p.y)), bool(p.smooth)))
        else:
            return "segment type %s" % p.type
        offs = []
    if offs:
        return "trailing off-curve points"
    return closed, complex(start.x, start.y), bool(start.smooth), out


def _contours_from_ufo_glyphs(glyphs):
    M = len(glyphs)
    if any(len(g.contours) != len(glyphs[0].contours) for g in glyphs):
        return "incompatible contour count"
    contours = []
    for ci in range(len(glyphs[0].contours)):
        per = [_ufo_contour(g.contours[ci]) for g in glyphs]
        for p in per:
            if isinstance(p, str):
                return p
        closed = per[0][0]
        count = len(per[0][3])
        if any(p[0] != closed or len(p[3]) != count for p in per):
            return "incompatible contour"
        if count == 0:
            continue  # an anchor
        for i in range(count):
            if any(p[3][i][0] != per[0][3][i][0] for p in per):
                return "incompatible segment types"
        segs = []
        flags = []
        cur = [p[1] for p in per]
        curflag = all(p[2] for p in per)
        for i in range(count):
            kind = per[0][3][i][0]
            if kind == 0:
                segs.append([(cur[m], per[m][3][i][1][0]) for m in range(M)])
            else:
                segs.append([(cur[m],) + per[m][3][i][1] for m in range(M)])
            flags.append(1 if curflag else 0)
            cur = [per[m][3][i][1][-1] for m in range(M)]
            curflag = all(p[3][i][2] for p in per)
        if closed:
            # The start point is the end of the last segment; its mark is the
            # one recorded there.
            flags[0] = 1 if all(p[3][-1][2] for p in per) else 0
        else:
            flags[0] = 0
        keep = [
            i
            for i, s in enumerate(segs)
            if not (len(s[0]) == 2 and all(p[0] == p[1] for p in s))
        ]
        if len(keep) != len(segs):
            segs = [segs[i] for i in keep]
            flags = [flags[i] for i in keep]
        contours.append((closed, segs, flags))
    return contours


def extract_designspace(path_or_doc, glyph_range=None):
    """Full masters of a designspace (or of a .glyphs file, through glyphsLib).

    Returns (arrays, master_locations (normalized), upem, skipped)."""
    from fontTools import varLib
    from fontTools.designspaceLib import DesignSpaceDocument

    if isinstance(path_or_doc, str) and path_or_doc.endswith(".glyphs"):
        import glyphsLib

        gsfont = glyphsLib.GSFont(path_or_doc)
        doc = glyphsLib.to_designspace(gsfont, minimize_glyphs_diffs=False)
        for s in doc.sources:
            if s.path is None:
                s.path = "/nonexistent/%s" % (s.filename or s.name)
    elif isinstance(path_or_doc, str):
        doc = DesignSpaceDocument.fromfile(path_or_doc)
    else:
        doc = path_or_doc
    data = varLib.load_designspace(doc, log_enabled=False)
    fonts = []
    locs = []
    for src, loc in zip(data.masters, data.normalized_master_locs):
        if getattr(src, "layerName", None):
            continue
        font = getattr(src, "font", None)
        if font is None:
            import ufoLib2

            font = ufoLib2.Font.open(src.path)
        fonts.append(font)
        locs.append(dict(loc))
    upem = fonts[0].info.unitsPerEm or 1000
    names = sorted(set.intersection(*[set(f.keys()) for f in fonts]))
    if glyph_range is not None:
        names = names[glyph_range[0] : glyph_range[1]]
    flat = _Flat(len(fonts))
    for name in names:
        glyphs = [f[name] for f in fonts]
        contours = _contours_from_ufo_glyphs(glyphs)
        if isinstance(contours, str):
            flat.skip(contours)
            continue
        if not contours:
            continue
        flat.add_glyph(name, contours)
    return flat.arrays(), locs, upem, flat.skipped


# ---------------------------------------------------------------------------
# conversion
# ---------------------------------------------------------------------------


def convert(arrays, tol, curves_to_quadratic=None):
    """Convert every cubic across masters. Adds n, qfirst, qlast to arrays."""
    if curves_to_quadratic is None:
        from fontTools.cu2qu.cu2qu import curves_to_quadratic
    pts = arrays["pts"]
    kind = arrays["kind"]
    S, M, _ = pts.shape
    n = np.zeros(S, dtype=np.int16)
    qfirst = np.full((S, M), complex(NAN, NAN), dtype=np.complex128)
    qlast = np.full((S, M), complex(NAN, NAN), dtype=np.complex128)
    tols = [tol] * M
    re = pts.real
    im = pts.imag
    for s in np.nonzero(kind == 1)[0]:
        curves = [
            [(re[s, m, k], im[s, m, k]) for k in range(4)] for m in range(M)
        ]
        splines = curves_to_quadratic(curves, tols)
        n[s] = len(splines[0]) - 2
        for m in range(M):
            sp = splines[m]
            qfirst[s, m] = complex(sp[1][0], sp[1][1])
            qlast[s, m] = complex(sp[-2][0], sp[-2][1])
    arrays["n"] = n
    arrays["qfirst"] = qfirst
    arrays["qlast"] = qlast
    return arrays
