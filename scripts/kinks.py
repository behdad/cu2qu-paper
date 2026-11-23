#!/usr/bin/env python3
"""How much cu2qu changes the angle at smooth joints of interpolated instances.

usage: kinks.py DATA_DIR [--jobs N] [--check N] [--only NAME] [--smooth-within DEGREES]

Run it with fontTools importable (PYTHONPATH=.../fonttools/Lib) and numpy.
The report goes to standard output, progress to standard error.

DATA_DIR holds, for every source, NAME.npz with its point-compatible cubic
masters, and NAME.npz.meta.json and NAME.npz.names.json beside it.  What is
read from them:

  pts[S, M, 4]  complex128  control points of S segments in M masters; a line
                            has p0, p1, nan, nan
  kind[S]       uint8       0 line, 1 cubic
  prev[S]       int64       the segment that ends where this one starts (the
                            contours are closed: the first segment of a
                            contour names the last)
  gid[S]        int32       index into the list of glyph names
  meta          "upem", "locations" (one per master), "skipped"

Nothing else in the files is read: the conversion is done here.

What is measured (the head of the report says it in full):

  * every cubic segment is converted across masters with
    fontTools.cu2qu.cu2qu.curves_to_quadratic, at 0.001 x units per em;
  * a joint is the start point of a segment, between the segment `prev` names
    (incoming) and the segment itself (outgoing);
  * the direction on a side of a joint is taken from the nearest control point
    of that side's segment that is not at the joint;
  * the kink is the angle between the two directions;
  * a joint is smooth if a side is a cubic and the kink is exactly 0 in every
    cubic master (--smooth-within makes that "at most DEGREES", for seeing how
    much the exactness matters; the default is exact);
  * for every pair of masters and t in 1/4, 1/2, 3/4, the instance
    (1 - t) x master i + t x master j is taken of the cubic masters, of the
    quadratic masters and of the quadratic masters with rounded off-curve
    points, giving kinks K_c, K_q and K_r at each smooth joint;
  * change by conversion = max over instances of |K_q - K_c|, change by
    rounding = max over instances of |K_r - K_q|;
  * besides, for what three instances per pair of masters can miss: the kink
    of each smooth joint in the quadratic masters themselves.
"""

import argparse
import glob
import hashlib
import json
import math
import multiprocessing
import os
import sys
import time

import numpy as np

TOL_PER_EM = 0.001
LAMBDAS = (0.25, 0.5, 0.75)
THRESHOLDS = (0.1, 0.5, 1.0)
# In an instance, a control point this close to the joint (font units) is at
# the joint.  The report says how far the nearest miss on either side is.
EPS = 1e-9
# Master coordinates are tested as whole multiples of 1/UNIT, which they are.
UNIT = 1000
NAN = float("nan")
CNAN = complex(NAN, NAN)
VERSIONS = "cqr"  # cubic, quadratic, quadratic with rounded off-curve points


# ---------------------------------------------------------------------------
# conversion
# ---------------------------------------------------------------------------


def _convert_chunk(job):
    """Convert cubics across masters; job = (coords[C, M, 4, 2], tolerance).

    Returns the number of quadratics of each (-1 if cu2qu raised), the
    off-curve points [sum n, M, 2] in order, the same rounded with otRound,
    what was raised, and how many results did not keep the end points.
    """
    from fontTools.cu2qu.cu2qu import curves_to_quadratic
    from fontTools.misc.roundTools import otRound

    coords, tol = job
    M = coords.shape[1]
    max_errors = [tol] * M
    ns = []
    flat = []
    flat_rounded = []
    raised = {}
    moved_ends = 0
    for curves in coords.tolist():
        try:
            splines = curves_to_quadratic(curves, max_errors)
        except Exception as e:  # counted and reported
            key = type(e).__name__
            raised[key] = raised.get(key, 0) + 1
            ns.append(-1)
            continue
        count = len(splines[0])
        for curve, spline in zip(curves, splines):
            if len(spline) != count:
                raise AssertionError("masters got different numbers of points")
            if (
                spline[0][0] != curve[0][0]
                or spline[0][1] != curve[0][1]
                or spline[-1][0] != curve[3][0]
                or spline[-1][1] != curve[3][1]
            ):
                moved_ends += 1
        ns.append(count - 2)
        for k in range(1, count - 1):
            for spline in splines:
                x, y = spline[k]
                flat.append(x)
                flat.append(y)
                flat_rounded.append(otRound(x))
                flat_rounded.append(otRound(y))
    return (
        np.array(ns, dtype=np.int64),
        np.array(flat, dtype=np.float64).reshape(-1, M, 2),
        np.array(flat_rounded, dtype=np.float64).reshape(-1, M, 2),
        raised,
        moved_ends,
    )


def convert(pts, cubics, tol, pool):
    """Convert the cubic segments with indices `cubics`, in that order."""
    M = pts.shape[1]
    size = max(200, 60000 // M)

    def jobs():
        for lo in range(0, len(cubics), size):
            p = pts[cubics[lo : lo + size]]
            yield np.stack([p.real, p.imag], axis=-1), tol

    results = pool.imap(_convert_chunk, jobs()) if pool else map(_convert_chunk, jobs())
    ns, qs, rs = [], [], []
    raised = {}
    moved_ends = 0
    for n, q, r, e, moved in results:
        ns.append(n)
        qs.append(q[..., 0] + 1j * q[..., 1])
        rs.append(r[..., 0] + 1j * r[..., 1])
        for k, v in e.items():
            raised[k] = raised.get(k, 0) + v
        moved_ends += moved
    if not ns:
        empty = np.zeros((0, M), dtype=complex)
        return np.zeros(0, dtype=np.int64), empty, empty, raised, moved_ends
    return np.concatenate(ns), np.concatenate(qs), np.concatenate(rs), raised, moved_ends


# ---------------------------------------------------------------------------
# the definitions, one number at a time (cross-check and worst-joint details)
# ---------------------------------------------------------------------------


def scalar_controls(pts, kind, seg):
    """Control points of a segment in every master: [M][2 or 4] complex."""
    k = 4 if kind[seg] == 1 else 2
    return [[complex(p) for p in pts[seg, m, :k]] for m in range(pts.shape[1])]


def scalar_quadratic(controls, tol, rounded):
    """The segment after conversion: [M][n + 2] complex (a line as it is)."""
    from fontTools.cu2qu.cu2qu import curves_to_quadratic
    from fontTools.misc.roundTools import otRound

    if len(controls[0]) == 2:
        return controls
    splines = curves_to_quadratic(
        [[(p.real, p.imag) for p in c] for c in controls], [tol] * len(controls)
    )
    out = []
    for c, spline in zip(controls, splines):
        off = [complex(x, y) for x, y in spline[1:-1]]
        if rounded:
            off = [complex(otRound(p.real), otRound(p.imag)) for p in off]
        out.append([c[0]] + off + [c[-1]])
    return out


def scalar_used(incoming, outgoing):
    """Which control point gives the direction on each side: two indices.

    Counted from the joint: 1 is the control point next to it.  None if all
    of the segment is at the joint.
    """
    joint = outgoing[0]
    used_in = used_out = None
    for k, p in enumerate(reversed(incoming[:-1])):
        if abs(p - joint) > EPS:
            used_in = k + 1
            break
    for k, p in enumerate(outgoing[1:]):
        if abs(p - joint) > EPS:
            used_out = k + 1
            break
    return used_in, used_out


def scalar_kink(incoming, outgoing):
    """Kink in degrees between a segment ending at a joint and one starting there."""
    joint = outgoing[0]
    used_in, used_out = scalar_used(incoming, outgoing)
    if used_in is None or used_out is None:
        return NAN
    din = joint - incoming[-1 - used_in]
    dout = outgoing[used_out] - joint
    w = din.conjugate() * dout
    return abs(math.degrees(math.atan2(w.imag, w.real)))


def scalar_joint(pts, kind, prev, s, instances, tol):
    """K_c, K_q, K_r of the joint at the start of segment s, per instance."""
    a = int(prev[s])
    sides = []
    for seg in (a, s):
        c = scalar_controls(pts, kind, seg)
        sides.append((c, scalar_quadratic(c, tol, False), scalar_quadratic(c, tol, True)))
    out = ([], [], [])
    for i, j, lam in instances:
        for which in range(3):
            A = [(1.0 - lam) * p + lam * q for p, q in zip(sides[0][which][i], sides[0][which][j])]
            B = [(1.0 - lam) * p + lam * q for p, q in zip(sides[1][which][i], sides[1][which][j])]
            out[which].append(scalar_kink(A, B))
    return out


# ---------------------------------------------------------------------------
# one source
# ---------------------------------------------------------------------------


def rowmax(x):
    """Row-wise maximum that ignores NaN: (max, argmax); NaN for all-NaN rows."""
    y = np.where(np.isnan(x), -np.inf, x)
    arg = y.argmax(axis=1)
    mx = y[np.arange(len(y)), arg]
    mx[np.isneginf(mx)] = NAN
    return mx, arg


def new_diag():
    return dict(nearest=math.inf, farthest_at_joint=0.0, undefined=0)


def fmt_points(ps):
    return " ".join("(%g, %g)" % (p.real, p.imag) for p in ps)


class Source:
    def __init__(self, path, smooth_within):
        self.name = os.path.basename(path)[: -len(".npz")]
        self.smooth_within = smooth_within
        data = np.load(path)
        self.pts = data["pts"]
        self.kind = data["kind"]
        self.prev = data["prev"].astype(np.int64)
        self.gid = data["gid"]
        with open(path + ".meta.json") as f:
            self.meta = json.load(f)
        with open(path + ".names.json") as f:
            self.names = json.load(f)
        self.S, self.M, four = self.pts.shape
        if four != 4:
            raise SystemExit("%s: unexpected shape of pts" % path)
        self.tol = TOL_PER_EM * self.meta["upem"]
        self.cub = self.kind == 1
        line = self.kind == 0
        if not np.all(self.cub | line):
            raise SystemExit("%s: segment kinds other than 0 and 1" % path)
        finite = np.isfinite(self.pts.real) & np.isfinite(self.pts.imag)
        if not (finite[:, :, :2].all() and finite[self.cub].all()):
            raise SystemExit("%s: missing control points" % path)
        if finite[line][:, :, 2:].any():
            raise SystemExit("%s: a line with more than two points" % path)
        # Exact arithmetic below relies on this.
        v = self.pts[finite]
        k = np.concatenate([v.real, v.imag]) * UNIT
        if np.abs(k - np.rint(k)).max() > 1e-6:
            raise SystemExit("%s: coordinates are not multiples of 1/%d" % (path, UNIT))
        if np.abs(k).max() > 2e7:
            raise SystemExit("%s: coordinates too large for the exact test" % path)
        frac = (v.real != np.rint(v.real)) | (v.imag != np.rint(v.imag))
        self.fractional_points = int(frac.sum())
        self.start = self.pts[:, :, 0]
        self.end = np.where(self.cub[:, None], self.pts[:, :, 3], self.pts[:, :, 1])
        pairs = [(i, j) for i in range(self.M) for j in range(i + 1, self.M)]
        self.instances = [(i, j, lam) for i, j in pairs for lam in LAMBDAS]
        self.II = np.array([i for i, j, lam in self.instances], dtype=np.int64)
        self.JJ = np.array([j for i, j, lam in self.instances], dtype=np.int64)
        self.LAM = np.array([lam for i, j, lam in self.instances], dtype=np.float64)
        # The masters themselves, in the form of instances: t = 0 of (m, m).
        self.masters = [(m, m, 0.0) for m in range(self.M)]
        self.master_instances = (np.arange(self.M), np.arange(self.M), np.zeros(self.M))

    # -- conversion ---------------------------------------------------------

    def run_conversion(self, pool):
        cub, M = self.cub, self.M
        cubics = np.nonzero(cub)[0]
        n, q, r, self.raised, self.moved_ends = convert(self.pts, cubics, self.tol, pool)
        self.n = np.zeros(self.S, dtype=np.int64)  # quadratics per segment; 0: line, -1: failed
        self.n[cubics] = n
        self.failed = self.n < 0
        count_q = np.where(self.failed, 0, self.n)
        off_q = np.concatenate([[0], np.cumsum(count_q)])
        if off_q[-1] != len(q):
            raise AssertionError("off-curve points do not add up")
        # Three versions of every segment, as (offset, count, table) of the
        # points strictly between its end points: the cubic's two control
        # points, the spline's off-curve points, and those rounded.
        count_c = np.where(cub, 2, 0)
        off_c = np.concatenate([[0], np.cumsum(count_c)])
        table_c = np.ascontiguousarray(self.pts[cub][:, :, 1:3].transpose(0, 2, 1)).reshape(-1, M)
        self.versions = {
            "c": (off_c, count_c, table_c),
            "q": (off_q, count_q, q),
            "r": (off_q, count_q, r),
        }
        self.single = cub & (self.n == 1)

    # -- smooth joints ------------------------------------------------------

    def find_smooth(self):
        """Classify the joints; a joint has the index of its outgoing segment."""
        pts, cub, prev, start, end = self.pts, self.cub, self.prev, self.start, self.end
        has_prev = prev >= 0
        pv = np.where(has_prev, prev, 0)
        joined = has_prev & np.all(start == end[pv], axis=1)
        self.joints = int(has_prev.sum())
        self.broken = int((has_prev & ~joined).sum())
        curved = np.nonzero(joined & (cub | cub[pv]))[0]
        self.curved_joints = len(curved)

        def units(z):
            return np.rint(z.real * UNIT) + 1j * np.rint(z.imag * UNIT)

        smooth_parts = []
        upright_parts = []
        undefined = 0
        undefined_otherwise_smooth = 0
        for lo in range(0, len(curved), 500000):
            js = curved[lo : lo + 500000]
            a = prev[js]
            J = units(start[js])
            # Incoming side: p2, else p1, else p0 of a cubic; p0 of a line.
            A = pts[a]
            ca = cub[a][:, None]
            din = J - units(np.where(ca, A[:, :, 2], A[:, :, 0]))
            din = np.where(din == 0, J - units(np.where(ca, A[:, :, 1], A[:, :, 0])), din)
            din = np.where(din == 0, J - units(A[:, :, 0]), din)
            # Outgoing side: p1, else p2, else p3 of a cubic; p1 of a line.
            B = pts[js]
            cs = cub[js][:, None]
            dout = units(B[:, :, 1]) - J
            dout = np.where(dout == 0, units(np.where(cs, B[:, :, 2], B[:, :, 1])) - J, dout)
            dout = np.where(dout == 0, units(np.where(cs, B[:, :, 3], B[:, :, 1])) - J, dout)
            # Whole numbers of at most 2e7 thousandths: products and their
            # sums are exact in double precision.
            cross = din.real * dout.imag - din.imag * dout.real
            dot = din.real * dout.real + din.imag * dout.imag
            defined = (din != 0) & (dout != 0)
            if self.smooth_within > 0:
                kink = np.degrees(np.abs(np.arctan2(cross, dot)))
                smooth_m = defined & (kink <= self.smooth_within)
            else:
                smooth_m = defined & (cross == 0) & (dot > 0)
            some_undefined = ~defined.all(axis=1)
            undefined += int(some_undefined.sum())
            undefined_otherwise_smooth += int((some_undefined & (smooth_m | ~defined).all(axis=1)).sum())
            is_smooth = smooth_m.all(axis=1)
            smooth_parts.append(js[is_smooth])
            # Both directions vertical in every master, or both horizontal in
            # every master.
            vertical = ((din.real == 0) & (dout.real == 0)).all(axis=1)
            horizontal = ((din.imag == 0) & (dout.imag == 0)).all(axis=1)
            upright_parts.append((vertical | horizontal)[is_smooth])
        smooth = np.concatenate(smooth_parts) if smooth_parts else np.zeros(0, dtype=np.int64)
        upright = np.concatenate(upright_parts) if upright_parts else np.zeros(0, dtype=bool)
        self.undefined_master = undefined
        self.undefined_master_otherwise_smooth = undefined_otherwise_smooth
        unconverted = self.failed[smooth] | self.failed[prev[smooth]]
        self.smooth_unconverted = int(unconverted.sum())
        self.smooth = smooth[~unconverted]
        self.upright = upright[~unconverted]

    # -- instances ----------------------------------------------------------

    def directions(self, version, seg, incoming, Ji, diag, instances):
        """Nearest control point not at the joint, minus the joint.

        For every joint (row; `seg` is its segment on this side) and instance
        (column).  NaN where the whole segment is at the joint.  Also returns
        where the control point next to the joint was at the joint.
        """
        off, count, table = self.versions[version]
        II, JJ, LAM = instances
        far = self.start if incoming else self.end
        n = count[seg]

        def candidate(sub, k):
            sg = seg[sub]
            nn = n[sub]
            C = np.full((len(sub), self.M), CNAN)
            inside = k < nn
            index = (off[sg] + nn - 1 - k) if incoming else (off[sg] + k)
            C[inside] = table[index[inside]]
            last = k == nn
            C[last] = far[sg[last]]
            return C, k <= nn

        C, _ = candidate(np.arange(len(seg)), 0)
        d = C[:, II] * (1.0 - LAM) + C[:, JJ] * LAM - Ji
        length = np.abs(d)
        at_joint = length <= EPS
        if (~at_joint).any():
            diag["nearest"] = min(diag["nearest"], float(length[~at_joint].min()))
        r, c = np.nonzero(at_joint)
        k = 0
        while len(r):
            diag["farthest_at_joint"] = max(diag["farthest_at_joint"], float(length[r, c].max()))
            k += 1
            C, exists = candidate(r, k)
            rows = np.arange(len(r))
            dd = C[rows, II[c]] * (1.0 - LAM[c]) + C[rows, JJ[c]] * LAM[c] - Ji[r, c]
            d[r, c] = dd
            ll = np.abs(dd)
            length[r, c] = ll
            diag["undefined"] += int((~exists).sum())
            used = exists & (ll > EPS)
            if used.any():
                diag["nearest"] = min(diag["nearest"], float(ll[used].min()))
            again = exists & (ll <= EPS)
            r, c = r[again], c[again]
        return d, at_joint

    def turns(self, js, diag, instances=None):
        """Signed kink in degrees, [joint, instance], for the three versions.

        Also, for each version, where the direction on a side did not come
        from the control point next to the joint.  `instances` is (i, j, t)
        as three arrays; the default is the instances of the report.
        """
        if instances is None:
            instances = (self.II, self.JJ, self.LAM)
        II, JJ, LAM = instances
        a = self.prev[js]
        Jm = self.start[js]
        Ji = Jm[:, II] * (1.0 - LAM) + Jm[:, JJ] * LAM
        turn = {}
        further = {}
        for version in VERSIONS:
            din, f_in = self.directions(version, a, True, Ji, diag[version], instances)
            dout, f_out = self.directions(version, js, False, Ji, diag[version], instances)
            # din is (control point - joint): the incoming direction is -din.
            w = np.conj(-din) * dout
            turn[version] = np.degrees(np.arctan2(w.imag, w.real))
            further[version] = f_in | f_out
        return turn, further

    def measure(self):
        js_all = self.smooth
        N = len(js_all)
        I = len(self.instances)
        self.diag = {v: new_diag() for v in VERSIONS}
        keys = (
            "conv",  # change by conversion
            "rnd",  # change by rounding
            "Kc_at",  # the kinks at the instance of the change by conversion
            "Kq_at",
            "Kr_at",
            "conv_regular",  # change by conversion over the regular instances
            "conv_signed",  # the changes for signed kinks
            "rnd_signed",
            "Kc_max",  # largest kink of the cubic instances
            "Kc_master",  # largest kink in a cubic master, in floating point
            "conv_master",  # largest |K_q - K_c| in a master, no instance taken
            "Kc_in_master",  # the two kinks in that master
            "Kq_in_master",
        )
        res = {k: np.full(N, NAN) for k in keys}
        res["conv_arg"] = np.zeros(N, dtype=np.int64)
        res["conv_master_arg"] = np.zeros(N, dtype=np.int64)
        scratch = {v: new_diag() for v in VERSIONS}
        # Instances in which, for the cubic or for the quadratic, a direction
        # did not come from the control point next to the joint.
        res["irregular"] = np.zeros(N, dtype=np.int64)
        for v in VERSIONS:
            res["further_" + v] = np.zeros(N, dtype=bool)
            res["undefined_" + v] = np.zeros(N, dtype=np.int64)
        wrap = lambda x: np.abs((x + 180.0) % 360.0 - 180.0)
        block = max(1, 3000000 // max(I, 1))
        with np.errstate(invalid="ignore"):
            for lo in range(0, N if I else 0, block):
                sl = slice(lo, lo + block)
                js = js_all[sl]
                turn, further = self.turns(js, self.diag)
                K = {v: np.abs(turn[v]) for v in VERSIONS}
                rows = np.arange(len(js))
                change = np.abs(K["q"] - K["c"])
                conv, arg = rowmax(change)
                res["conv"][sl] = conv
                res["conv_arg"][sl] = arg
                res["Kc_at"][sl] = K["c"][rows, arg]
                res["Kq_at"][sl] = K["q"][rows, arg]
                res["Kr_at"][sl] = K["r"][rows, arg]
                res["rnd"][sl] = rowmax(np.abs(K["r"] - K["q"]))[0]
                res["Kc_max"][sl] = rowmax(K["c"])[0]
                irregular = further["c"] | further["q"]
                res["irregular"][sl] = irregular.sum(axis=1)
                res["conv_regular"][sl] = rowmax(np.where(irregular, NAN, change))[0]
                # With signed turning angles (not what the tables report): a
                # difference of angles, wrapped to +-180.
                res["conv_signed"][sl] = rowmax(wrap(turn["q"] - turn["c"]))[0]
                res["rnd_signed"][sl] = rowmax(wrap(turn["r"] - turn["q"]))[0]
                for v in VERSIONS:
                    res["further_" + v][sl] = further[v].any(axis=1)
                    res["undefined_" + v][sl] = np.isnan(turn[v]).sum(axis=1)
                # The masters themselves, with no interpolation.
                turn, _ = self.turns(js, scratch, self.master_instances)
                Kc, Kq = np.abs(turn["c"]), np.abs(turn["q"])
                res["Kc_master"][sl] = rowmax(Kc)[0]
                res["conv_master"][sl], arg = rowmax(np.abs(Kq - Kc))
                res["conv_master_arg"][sl] = arg
                res["Kc_in_master"][sl] = Kc[rows, arg]
                res["Kq_in_master"][sl] = Kq[rows, arg]
        res["single"] = self.single[js_all] | self.single[self.prev[js_all]]
        res["upright"] = self.upright
        J = self.start[js_all]
        res["fractional"] = ((J.real != np.rint(J.real)) | (J.imag != np.rint(J.imag))).any(axis=1)
        self.res = res

    # -- cross-check and details --------------------------------------------

    def cross_check(self, count):
        """Compare the array code with the scalar code on a sample of joints."""
        N = len(self.smooth)
        if not N or not self.instances or count <= 0:
            return 0, 0.0, 0
        res = self.res
        rng = np.random.default_rng(20261123)
        single = np.nonzero(res["single"])[0]
        other = np.nonzero(~res["single"])[0]
        pick = [
            rng.choice(g, size=min(len(g), (count + 1) // 2), replace=False) for g in (single, other)
        ]
        # ... and the joints that matter most or are the least ordinary.
        conv = np.where(np.isnan(res["conv"]), -1.0, res["conv"])
        rnd = np.where(np.isnan(res["rnd"]), -1.0, res["rnd"])
        odd = np.nonzero(res["further_c"] | res["further_q"] | res["further_r"])[0]
        pick += [np.argsort(-conv)[:20], np.argsort(-rnd)[:20], odd[:100]]
        master = np.where(np.isnan(res["conv_master"]), -1.0, res["conv_master"])
        pick.append(np.argsort(-master)[:20])
        sample = np.unique(np.concatenate(pick))
        scratch = {v: new_diag() for v in VERSIONS}
        with np.errstate(invalid="ignore"):
            # In the instances of the report, and in the masters themselves.
            cases = [
                (self.instances, self.turns(self.smooth[sample], scratch)[0]),
                (self.masters, self.turns(self.smooth[sample], scratch, self.master_instances)[0]),
            ]
        worst = 0.0
        nan_mismatch = 0
        for row, k in enumerate(sample):
            for instances, turn in cases:
                scalar = scalar_joint(self.pts, self.kind, self.prev, int(self.smooth[k]), instances, self.tol)
                for v, sc in zip(VERSIONS, scalar):
                    vec = np.abs(turn[v][row])
                    sc = np.array(sc)
                    nv, ns = np.isnan(vec), np.isnan(sc)
                    nan_mismatch += int((nv != ns).sum())
                    both = ~nv & ~ns
                    if both.any():
                        worst = max(worst, float(np.abs(vec[both] - sc[both]).max()))
        return len(sample), worst, nan_mismatch

    def joint_record(self, k):
        """Summary of smooth joint number k (index into self.smooth)."""
        s = int(self.smooth[k])
        a = int(self.prev[s])
        res = self.res
        i, j, lam = self.instances[int(res["conv_arg"][k])]
        return dict(
            source=self.name,
            glyph=self.names[int(self.gid[s])],
            s=s,
            a=a,
            n_in=int(self.n[a]),
            n_out=int(self.n[s]),
            i=i,
            j=j,
            lam=lam,
            conv=float(res["conv"][k]),
            rnd=float(res["rnd"][k]),
            Kc=float(res["Kc_at"][k]),
            Kq=float(res["Kq_at"][k]),
            Kr=float(res["Kr_at"][k]),
            irregular=int(res["irregular"][k]),
            conv_master=float(res["conv_master"][k]),
            master=int(res["conv_master_arg"][k]),
            Kc_in_master=float(res["Kc_in_master"][k]),
            Kq_in_master=float(res["Kq_in_master"][k]),
            single=bool(res["single"][k]),
        )

    def describe(self, rec):
        """A worst-joint paragraph: control points and kinks, recomputed."""
        s, a, i, j, lam = rec["s"], rec["a"], rec["i"], rec["j"], rec["lam"]
        lines = []
        add = lines.append
        add("source %s, glyph %s" % (self.name, rec["glyph"]))
        add(
            "joint at the start of segment %d (incoming segment %d); masters %d and %d, t = %g"
            % (s, a, i, j, lam)
        )
        locs = self.meta.get("locations") or []
        for m in (i, j):
            if m < len(locs):
                add("  master %d is at %s" % (m, json.dumps(locs[m], sort_keys=True)))
        ctrl = {seg: scalar_controls(self.pts, self.kind, seg) for seg in (a, s)}
        quad = {seg: scalar_quadratic(ctrl[seg], self.tol, False) for seg in (a, s)}
        rquad = {seg: scalar_quadratic(ctrl[seg], self.tol, True) for seg in (a, s)}
        for seg, side in ((a, "incoming"), (s, "outgoing")):
            if self.kind[seg] == 0:
                what = "line"
            else:
                what = "cubic, converted to %d quadratic%s" % (self.n[seg], "" if self.n[seg] == 1 else "s")
            add("  %s segment (%s):" % (side, what))
            for m in (i, j):
                add("    master %-2d  %s" % (m, fmt_points(ctrl[seg][m])))
                if self.kind[seg] == 1:
                    add("      quadratic  %s" % fmt_points(quad[seg][m]))
                    add("      rounded    %s" % fmt_points(rquad[seg][m]))
        inst = lambda P: [(1.0 - lam) * p + lam * q for p, q in zip(P[i], P[j])]
        add("  the instance:")
        for seg, side in ((a, "incoming"), (s, "outgoing")):
            add("    %s  %s" % (side, fmt_points(inst(ctrl[seg]))))
            if self.kind[seg] == 1:
                add("      quadratic  %s" % fmt_points(inst(quad[seg])))
                add("      rounded    %s" % fmt_points(inst(rquad[seg])))
        Kc = scalar_kink(inst(ctrl[a]), inst(ctrl[s]))
        Kq = scalar_kink(inst(quad[a]), inst(quad[s]))
        Kr = scalar_kink(inst(rquad[a]), inst(rquad[s]))
        add("  K_c = %.6f, K_q = %.6f, K_r = %.6f degrees (recomputed one number at a time)" % (Kc, Kq, Kr))
        add(
            "  the joint's change by conversion is %.6f; its change by rounding, over all instances, %.6f"
            % (rec["conv"], rec["rnd"])
        )
        used = [
            scalar_used(inst(ctrl[a]), inst(ctrl[s])),
            scalar_used(inst(quad[a]), inst(quad[s])),
        ]
        if any(u != 1 for pair in used for u in pair):
            add(
                "  in this instance a control point next to the joint is at the joint; the directions come from"
            )
            add(
                "    control point number (counting from the joint) %s/%s of the incoming/outgoing cubic, %s/%s"
                % (used[0] + used[1])
            )
            add("    of the quadratics")
        mc = ["%.4g" % scalar_kink(ctrl[a][m], ctrl[s][m]) for m in range(self.M)]
        add("  kink in each cubic master: %s" % " ".join(mc))
        mq = ["%.4g" % scalar_kink(quad[a][m], quad[s][m]) for m in range(self.M)]
        add("  kink in each quadratic master, not rounded: %s" % " ".join(mq))
        return "\n".join(lines)


def over(x, mask):
    """Counts above the thresholds and the largest, over x[mask], ignoring NaN."""
    v = x[mask]
    v = v[~np.isnan(v)]
    return [int((v > t).sum()) for t in THRESHOLDS] + [float(v.max()) if len(v) else NAN]


def process(path, pool, args):
    t0 = time.time()
    src = Source(path, args.smooth_within)
    src.run_conversion(pool)
    t1 = time.time()
    src.find_smooth()
    src.measure()
    t2 = time.time()
    checked = src.cross_check(args.check)
    res = src.res
    single = res["single"]
    other = ~single
    irregular = res["irregular"] > 0
    order = np.argsort(-np.where(np.isnan(res["conv"]), -1.0, res["conv"]), kind="stable")
    top_single = [src.joint_record(int(k)) for k in order[single[order]][:10]]
    top_other = [src.joint_record(int(k)) for k in order[other[order]][:10]]
    in_master = np.where(np.isnan(res["conv_master"]), -1.0, res["conv_master"])
    top_master = [src.joint_record(int(k)) for k in np.argsort(-in_master, kind="stable")[:10] if in_master[k] > THRESHOLDS[0]]
    out = dict(
        name=src.name,
        masters=src.M,
        instances=len(src.instances),
        glyphs=int(len(np.unique(src.gid))),
        names=len(src.names),
        segments=src.S,
        cubics=int(src.cub.sum()),
        singles=int(src.single.sum()),
        raised=src.raised,
        moved_ends=src.moved_ends,
        n_max=int(src.n.max()) if src.S else 0,
        joints=src.joints,
        broken=src.broken,
        curved_joints=src.curved_joints,
        undefined_master=src.undefined_master,
        undefined_master_otherwise_smooth=src.undefined_master_otherwise_smooth,
        smooth_unconverted=src.smooth_unconverted,
        smooth=len(src.smooth),
        smooth_single=int(single.sum()),
        smooth_single_both=int((src.single[src.smooth] & src.single[src.prev[src.smooth]]).sum()),
        smooth_upright=int(res["upright"].sum()),
        smooth_single_upright=int((res["upright"] & single).sum()),
        # counts above the thresholds and the largest, by population
        conv_single=over(res["conv"], single),
        rnd_single=over(res["rnd"], single),
        conv_other=over(res["conv"], other),
        rnd_other=over(res["rnd"], other),
        conv_regular_single=over(res["conv_regular"], single),
        conv_regular_other=over(res["conv_regular"], other),
        conv_signed_single=over(res["conv_signed"], single),
        rnd_signed_single=over(res["rnd_signed"], single),
        base_single=over(res["Kc_max"], single),
        base_other=over(res["Kc_max"], other),
        master_single=over(res["conv_master"], single),
        master_other=over(res["conv_master"], other),
        master_cubic=over(res["Kc_master"], np.ones(len(single), dtype=bool)),
        top_master=top_master,
        slanted_values=(res["conv"][single & ~res["upright"]], res["rnd"][single & ~res["upright"]]),
        open_starts=int((src.prev < 0).sum()),
        conv_single_slanted=over(res["conv"], single & ~res["upright"]),
        rnd_single_slanted=over(res["rnd"], single & ~res["upright"]),
        conv_single_upright=over(res["conv"], single & res["upright"]),
        rnd_single_upright=over(res["rnd"], single & res["upright"]),
        conv_upright=over(res["conv"], res["upright"]),
        rnd_upright=over(res["rnd"], res["upright"]),
        base_upright=over(res["Kc_max"], res["upright"]),
        rnd_single_fractional=over(res["rnd"], single & res["fractional"]),
        rnd_other_fractional=over(res["rnd"], other & res["fractional"]),
        fractional_single=int((single & res["fractional"]).sum()),
        fractional_other=int((other & res["fractional"]).sum()),
        irregular_single=int((irregular & single).sum()),
        irregular_other=int((irregular & other).sum()),
        irregular_pairs_single=int(res["irregular"][single].sum()),
        irregular_pairs_other=int(res["irregular"][other].sum()),
        nan_conv=int(np.isnan(res["conv"]).sum()),
        nan_rnd=int(np.isnan(res["rnd"]).sum()),
        undefined={v: int(res["undefined_" + v].sum()) for v in VERSIONS},
        further={v: int(res["further_" + v].sum()) for v in VERSIONS},
        further_single={v: int((res["further_" + v] & single).sum()) for v in VERSIONS},
        diag=src.diag,
        checked=checked,
        fractional_points=src.fractional_points,
        skipped=src.meta.get("skipped") or {},
        upem=src.meta["upem"],
        tol=src.tol,
        top_single=top_single,
        top_other=top_other,
        detail_single=src.describe(top_single[0]) if top_single else "",
        detail_other=src.describe(top_other[0]) if top_other else "",
        seconds=(t1 - t0, t2 - t1, time.time() - t2),
    )
    return out


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------


def deg(x):
    if x != x:
        return "-"
    if x == 0:
        return "0"
    if x < 0.001:
        return "%.1e" % x
    return "%.4f" % x


def table(header, rows, out):
    widths = [max(len(str(r[c])) for r in [header] + rows) for c in range(len(header))]
    for r in [header] + rows:
        cells = [str(r[0]).ljust(widths[0])] + [str(x).rjust(w) for x, w in zip(r[1:], widths[1:])]
        out("  " + "  ".join(cells).rstrip())
        if r is header:
            out("  " + "  ".join("-" * w for w in widths))


def merge_over(parts):
    """Combine results of over(): sum the counts, take the largest maximum."""
    counts = [sum(p[k] for p in parts) for k in range(len(THRESHOLDS))]
    maxima = [p[-1] for p in parts if p[-1] == p[-1]]
    return counts + [max(maxima) if maxima else NAN]


def pct(a, b):
    return "%.2f%%" % (100.0 * a / b) if b else "-"


def conv_key(rec):
    return -(rec["conv"] if rec["conv"] == rec["conv"] else -1.0)


def fonttools_identity():
    import fontTools
    from fontTools.cu2qu import cu2qu

    with open(cu2qu.__file__, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    branch = ""
    head = os.path.join(os.path.dirname(fontTools.__file__), "..", "..", ".git", "HEAD")
    try:
        with open(head) as f:
            branch = f.read().strip()
    except OSError:
        pass
    return fontTools.version, cu2qu.__file__, cu2qu.COMPILED, digest, branch


def report(results, args, seconds, out):
    version, cu2qu_file, compiled, digest, branch = fonttools_identity()
    tot = lambda k: sum(r[k] for r in results)
    merged = lambda k: merge_over([r[k] for r in results])
    exact = not args.smooth_within > 0

    out("Kinks at smooth joints of interpolated instances: cubic masters, cu2qu, rounding")
    out("=" * 82)
    out("")
    out("data       %s (%d sources)" % (os.path.abspath(args.data_dir), len(results)))
    out("fontTools  %s, %s" % (version, cu2qu_file))
    out("           compiled: %s; sha256 of cu2qu.py %s" % (compiled, digest))
    if branch:
        out("           HEAD of that checkout: %s" % branch)
    out("python     %s, numpy %s" % (sys.version.split()[0], np.__version__))
    out("")
    out("What is measured")
    out("")
    if exact:
        smooth_text = [
            "smooth       a side is a cubic, and in every master both directions are defined, their",
            "             cross product is exactly 0 and their dot product is positive. Every coordinate",
            "             in the data is a multiple of 1/%d and the test is done in whole 1/%d units:" % (UNIT, UNIT),
            "             it is exact, with no tolerance, also where coordinates are not integers.",
        ]
    else:
        smooth_text = [
            "smooth       NOT THE EXACT TEST: a side is a cubic, and in every master both directions",
            "             are defined and the kink is at most %g degrees (--smooth-within)." % args.smooth_within,
        ]
    for text in (
        [
            "tolerance    %g x units per em, the same for all masters" % TOL_PER_EM,
            "conversion   curves_to_quadratic([the segment in each master], [tolerance] * M) for every",
            "             cubic segment; n = points per master - 2 (n = 1: a single quadratic). Lines stay.",
            "joint        the start point of a segment, between the segment before it in its contour",
            "             (incoming) and the segment itself (outgoing). "
            + ("All contours in the data are closed." if not tot("open_starts") else "%d segments start an open contour: no joint there." % tot("open_starts")),
            "direction    incoming: the joint minus the nearest earlier control point of the incoming",
            "             segment that is not at the joint (cubic: p2, else p1, else p0; quadratic spline:",
            "             the last off-curve point, else the ones before it in turn, else p0; line: p0).",
            "             Outgoing: likewise forwards (p1, else p2, else p3; first off-curve point, ...).",
            "             Undefined if the whole segment is at the joint.",
            "             'At the joint': in masters, equal; in instances, within %g units." % EPS,
            "kink K       the angle between the two directions, 0 to 180 degrees, without a sign",
        ]
        + smooth_text
        + [
            "instances    (1 - t) x master i + t x master j for every pair i < j and t in 1/4, 1/2, 3/4;",
            "             all pairs, no sampling",
            "K_c K_q K_r  the kink in the instance of the cubic masters; of the quadratic masters; of the",
            "             quadratic masters with every off-curve coordinate rounded by otRound (on-curve",
            "             points are the cubic's and are not rounded)",
            "changes      by conversion: max over instances of |K_q - K_c|",
            "             by rounding:   max over instances of |K_r - K_q|",
            "             (an instance where one of the two kinks is undefined is left out)",
        ]
    ):
        out("  " + text)
    out("")

    # 1
    out("1. Sources and conversion")
    out("")
    header = ["source", "masters", "instances", "glyphs", "segments", "cubics", "single (n=1)", "share"]
    rows = [
        [r["name"], r["masters"], r["instances"], r["glyphs"], r["segments"], r["cubics"], r["singles"], pct(r["singles"], r["cubics"])]
        for r in results
    ]
    rows.append(["total", "", "", tot("glyphs"), tot("segments"), tot("cubics"), tot("singles"), pct(tot("singles"), tot("cubics"))])
    table(header, rows, out)
    out("")
    out("  glyphs: those with at least one segment in the data. instances: pairs of masters x 3.")
    out("  single: cubic segments that became a single quadratic (in every master: n is shared).")
    out("")

    # 2
    out("2. Joints")
    out("")
    header = ["source", "joints", "with a cubic side", "smooth", "not H/V", "next to a single", "share", "not H/V", "undefined"]

    def row2(name, g):
        return [
            name,
            g("joints"),
            g("curved_joints"),
            g("smooth"),
            g("smooth") - g("smooth_upright"),
            g("smooth_single"),
            pct(g("smooth_single"), g("smooth")),
            g("smooth_single") - g("smooth_single_upright"),
            g("undefined_master"),
        ]

    rows = [row2(r["name"], r.__getitem__) for r in results] + [row2("total", tot)]
    table(header, rows, out)
    out("")
    out("  next to a single: smooth joints with a single quadratic on at least one side (on both")
    out("    sides: %d of the %d)." % (tot("smooth_single_both"), tot("smooth_single")))
    cu, ru, bu = merged("conv_upright"), merged("rnd_upright"), merged("base_upright")
    out("  not H/V: of the column before, the joints that are not H/V. A joint is H/V if both its")
    out("    directions are horizontal in every master, or both vertical in every master: %d of the" % tot("smooth_upright"))
    out("    %d smooth joints. In an instance the control points next to such a joint stay on that" % tot("smooth"))
    out("    line through it, converted or not, so its kink is 0, or 180 where the two sides fold onto")
    out("    each other, unless a control point reaches the joint or the joint is not on integers.")
    out("    Of the H/V joints, %d have a K_c above 0.1 degree in some instance, %d a change by" % (bu[0], cu[0]))
    out("    conversion above 0.1, and %d a change by rounding above 0.1." % ru[0])
    out("  undefined: joints with a cubic side that are not counted as smooth because a direction")
    out("    is undefined in some master (a side has zero length there); %d of the %d are" % (tot("undefined_master_otherwise_smooth"), tot("undefined_master")))
    out("    smooth in all the masters where both directions are defined.")
    out("")

    # 3
    out("3. Smooth joints next to a single quadratic: joints whose change is above 0.1, 0.5, 1 degree")
    out("")
    header = ["source", "joints", "conv >0.1", ">0.5", ">1", "largest", "round >0.1", ">0.5", ">1", "largest"]
    rows = []
    for r in results:
        c, q = r["conv_single"], r["rnd_single"]
        rows.append([r["name"], r["smooth_single"]] + c[:3] + [deg(c[3])] + q[:3] + [deg(q[3])])
    c, q = merged("conv_single"), merged("rnd_single")
    rows.append(["total", tot("smooth_single")] + c[:3] + [deg(c[3])] + q[:3] + [deg(q[3])])
    table(header, rows, out)
    out("")
    n = tot("smooth_single")
    out("  total, as shares of the %d joints: conversion %s, %s, %s; rounding %s, %s, %s" % ((n,) + tuple(pct(x, n) for x in c[:3] + q[:3])))
    slanted = n - tot("smooth_single_upright")
    c2, q2 = merged("conv_single_slanted"), merged("rnd_single_slanted")
    out("  of the %d that are not H/V (see 2): conversion %d, %d, %d (%s, %s, %s), largest %s;" % ((slanted,) + tuple(c2[:3]) + tuple(pct(x, slanted) for x in c2[:3]) + (deg(c2[3]),)))
    out("    rounding %d, %d, %d (%s, %s, %s), largest %s" % (tuple(q2[:3]) + tuple(pct(x, slanted) for x in q2[:3]) + (deg(q2[3]),)))
    values = [np.concatenate([r["slanted_values"][k] for r in results]) for k in (0, 1)]
    values = [v[~np.isnan(v)] for v in values]
    if len(values[0]) and len(values[1]):
        ranks = (50, 90, 99)
        pc = [[deg(float(x)) for x in np.percentile(v, ranks)] for v in values]
        out("    50th, 90th, 99th percentile of the change over these joints: conversion %s, %s, %s;" % tuple(pc[0]))
        out("    rounding %s, %s, %s degrees" % tuple(pc[1]))
    c2, q2 = merged("conv_single_upright"), merged("rnd_single_upright")
    out("  of the %d that are H/V: conversion %d, %d, %d, largest %s; rounding %d, %d, %d, largest %s" % ((n - slanted,) + tuple(c2[:3]) + (deg(c2[3]),) + tuple(q2[:3]) + (deg(q2[3]),)))
    cr = merged("conv_regular_single")
    if tot("irregular_single"):
        out("  %d of these joints have an instance (%d joint-instance pairs) in which a control point next" % (tot("irregular_single"), tot("irregular_pairs_single")))
        out("    to the joint is at the joint, for the cubic or the quadratic, so that a direction comes")
        out("    from a point further away. Leaving those instances out, the change by conversion is")
        out("    above 0.1, 0.5, 1 degree at %d, %d, %d joints, largest %s." % (cr[0], cr[1], cr[2], deg(cr[3])))
    else:
        out("  none of these joints has an instance in which a control point next to the joint is at the")
        out("    joint, for the cubic or the quadratic (compare 4).")
    cs, qs = merged("conv_signed_single"), merged("rnd_signed_single")
    out("  for comparison only: if the kinks carried a sign (a left or a right turn) and the changes")
    out("    were differences of signed angles, the totals would be: conversion %d, %d, %d, largest %s;" % (cs[0], cs[1], cs[2], deg(cs[3])))
    out("    rounding %d, %d, %d, largest %s." % (qs[0], qs[1], qs[2], deg(qs[3])))
    b = merged("base_single")
    out("  for scale, the cubic instances themselves: K_c is above 0.1, 0.5, 1 degree in some instance")
    out("    at %d, %d, %d of these joints, largest %s." % (b[0], b[1], b[2], deg(b[3])))
    b = merged("master_single")
    out("  not an instance, but what the instances are made of: in some master itself, |K_q - K_c| (not")
    out("    rounded) is above 0.1, 0.5, 1 degree at %d, %d, %d of these joints, largest %s (see 5)." % (b[0], b[1], b[2], deg(b[3])))
    out("")

    # 4
    out("4. Smooth joints not next to a single quadratic")
    out("")
    header = ["source", "joints", "largest conv", "regular instances", "irregular joints", "round >0.1", ">0.5", ">1", "largest"]
    rows = []
    for r in results:
        q = r["rnd_other"]
        rows.append(
            [r["name"], r["smooth"] - r["smooth_single"], deg(r["conv_other"][3]), deg(r["conv_regular_other"][3]), r["irregular_other"]]
            + q[:3]
            + [deg(q[3])]
        )
    co, cr, q = merged("conv_other"), merged("conv_regular_other"), merged("rnd_other")
    rows.append(["total", tot("smooth") - tot("smooth_single"), deg(co[3]), deg(cr[3]), tot("irregular_other")] + q[:3] + [deg(q[3])])
    table(header, rows, out)
    out("")
    out("  largest conv: the largest change by conversion, in degrees, over all instances.")
    out("  regular instances: the same, leaving out the instances in which a control point next to")
    out("    the joint is at the joint (the handles of the two masters cancel), for the cubic or the")
    out("    quadratic, so that a direction comes from a point further away.")
    out("  irregular joints: the joints that have such an instance (%d joint-instance pairs in all)." % tot("irregular_pairs_other"))
    out("  joints with a change by conversion above 0.1, 0.5, 1 degree: %d, %d, %d over all instances;" % (co[0], co[1], co[2]))
    out("    %d, %d, %d over regular instances." % (cr[0], cr[1], cr[2]))
    b = merged("base_other")
    out("  for scale, the cubic instances themselves: K_c is above 0.1, 0.5, 1 degree in some instance")
    out("    at %d, %d, %d of these joints, largest %s." % (b[0], b[1], b[2], deg(b[3])))
    b = merged("master_other")
    out("  in some master itself, |K_q - K_c| (not rounded) is above 0.1, 0.5, 1 degree at %d, %d, %d" % (b[0], b[1], b[2]))
    out("    of these joints, largest %s (see 5)." % deg(b[3]))
    out("")

    # 5
    out("5. The worst joints for the change by conversion")
    out("")
    header = ["source", "glyph", "segment", "n in/out", "masters", "t", "K_c", "K_q", "K_r", "conversion", "rounding", "irregular"]

    def row5(rec):
        return [
            rec["source"],
            rec["glyph"],
            rec["s"],
            "%d/%d" % (rec["n_in"], rec["n_out"]),
            "%d,%d" % (rec["i"], rec["j"]),
            "%g" % rec["lam"],
            deg(rec["Kc"]),
            deg(rec["Kq"]),
            deg(rec["Kr"]),
            deg(rec["conv"]),
            deg(rec["rnd"]),
            rec["irregular"],
        ]

    for key, detail, title in (
        ("top_single", "detail_single", "Next to a single quadratic"),
        ("top_other", "detail_other", "Not next to a single quadratic"),
    ):
        ranked = sorted((rec for r in results for rec in r[key]), key=conv_key)
        if not ranked:
            continue
        out("  %s, the worst:" % title)
        out("")
        for r in results:
            if r["name"] == ranked[0]["source"]:
                for text in r[detail].split("\n"):
                    out("    " + text)
        out("")
        out("  %s, the ten worst:" % title)
        out("")
        table(header, [row5(rec) for rec in ranked[:10]], out)
        out("")
    out("  The worst of each source, of all its smooth joints:")
    out("")
    rows = []
    for r in results:
        best = sorted(r["top_single"][:1] + r["top_other"][:1], key=conv_key)
        if best:
            rows.append(row5(best[0]))
    table(header, rows, out)
    out("")
    out("  segment: index of the outgoing segment in the source's arrays. n in/out: quadratics of the")
    out("  incoming and outgoing segment (0: a line). K_c, K_q, K_r: at the instance of the largest")
    out("  |K_q - K_c|. rounding: the joint's change by rounding, over all instances. irregular: the")
    out("  joint's instances in which a control point next to the joint is at the joint.")
    out("")
    ranked = sorted((rec for r in results for rec in r["top_master"]), key=lambda rec: -rec["conv_master"])
    out("  Smooth joints where the conversion changes the kink of a master itself, with no instance")
    out("  taken, by more than 0.1 degree (not rounded). At most ten, the worst first:")
    out("")
    if ranked:
        header = ["source", "glyph", "segment", "n in/out", "master", "K_c there", "K_q there", "conversion", "next to a single"]
        rows = [
            [
                rec["source"],
                rec["glyph"],
                rec["s"],
                "%d/%d" % (rec["n_in"], rec["n_out"]),
                rec["master"],
                deg(rec["Kc_in_master"]),
                deg(rec["Kq_in_master"]),
                deg(rec["conv"]),
                "yes" if rec["single"] else "no",
            ]
            for rec in ranked[:10]
        ]
        table(header, rows, out)
        out("")
        out("  K_c there, K_q there: the joint's kink in that cubic master and in that quadratic master.")
        out("  conversion: the joint's change by conversion over the instances, as in the tables above.")
    else:
        out("  none")
    out("")

    # notes
    out("Notes")
    out("")
    raised = {}
    for r in results:
        for k, v in r["raised"].items():
            raised[k] = raised.get(k, 0) + v
    out("  conversions that raised an exception: %s" % (json.dumps(raised, sort_keys=True) if raised else "none"))
    out("  conversions that moved an end point: %d; smooth joints left out for a failed conversion: %d" % (tot("moved_ends"), tot("smooth_unconverted")))
    out("  largest n: %d" % max(r["n_max"] for r in results))
    out("  joints whose two segments do not meet exactly in every master: %d" % tot("broken"))
    out("  largest kink of a smooth joint in a cubic master, as floating point computes it: %s degrees" % deg(merged("master_cubic")[3]))
    pairs = sum(r["smooth"] * r["instances"] for r in results)
    out("  joint-instance pairs measured: %d" % pairs)
    und = {v: sum(r["undefined"][v] for r in results) for v in VERSIONS}
    out("    with an undefined kink: cubic %d, quadratic %d, rounded %d" % (und["c"], und["q"], und["r"]))
    out("  smooth joints with no instance where the change is defined: conversion %d, rounding %d" % (tot("nan_conv"), tot("nan_rnd")))
    fur = {v: sum(r["further"][v] for r in results) for v in VERSIONS}
    furs = {v: sum(r["further_single"][v] for r in results) for v in VERSIONS}
    out("  smooth joints where, in some instance, a direction came from a control point further away")
    out("    because the one next to the joint was at the joint: cubic %d, quadratic %d, rounded %d" % (fur["c"], fur["q"], fur["r"]))
    out("    (of which next to a single quadratic: %d, %d, %d)" % (furs["c"], furs["q"], furs["r"]))
    out("  the %g threshold for 'at the joint' in instances: over all directions taken," % EPS)
    for v, label in (("c", "cubic"), ("q", "quadratic"), ("r", "rounded")):
        near = min(r["diag"][v]["nearest"] for r in results)
        far = max(r["diag"][v]["farthest_at_joint"] for r in results)
        out("    %-9s shortest distance counted as away from the joint %.3g, longest counted as at it %.3g" % (label, near, far))
    frac = [r["name"] for r in results if r["fractional_points"]]
    out("  control points with a coordinate that is not an integer: %d, in %d sources" % (tot("fractional_points"), len(frac)))
    fs, fo = merged("rnd_single_fractional"), merged("rnd_other_fractional")
    out("    smooth joints at such a point in some master: %d next to a single quadratic (change by" % tot("fractional_single"))
    out("    rounding above 0.1 degree at %d of them), %d others (at %d). There, rounding moves the" % (fs[0], tot("fractional_other"), fo[0]))
    out("    off-curve points onto integers while the joint stays off them.")
    skipped = ["%s %s" % (r["name"], json.dumps(r["skipped"], sort_keys=True)) for r in results if r["skipped"]]
    out("  glyphs skipped when the masters were extracted (from the meta files): %s" % ("; ".join(skipped) if skipped else "none"))
    extra = ["%s %d" % (r["name"], r["names"] - r["glyphs"]) for r in results if r["names"] != r["glyphs"]]
    out("  glyph names without segments: %s" % ("; ".join(extra) if extra else "none"))
    out("  units per em: %s; tolerance: %s" % (sorted(set(r["upem"] for r in results)), sorted(set(r["tol"] for r in results))))
    nchecked = sum(r["checked"][0] for r in results)
    worst_diff = max([r["checked"][1] for r in results] or [0.0])
    mismatch = sum(r["checked"][2] for r in results)
    out("  cross-check: K_c, K_q, K_r of %d joints, in every instance and every master, recomputed one" % nchecked)
    out("    number at a time from the definitions, converting again: largest difference %.3g degrees," % worst_diff)
    out("    %d kinks defined in one computation and not in the other%s" % (mismatch, "" if worst_diff < 1e-6 and not mismatch else "  *** CROSS-CHECK FAILED ***"))
    out("")
    out(
        "run time: %.0f seconds with %d processes (conversion %.0f, joints and instances %.0f, cross-check and details %.0f)"
        % (
            seconds,
            args.jobs,
            sum(r["seconds"][0] for r in results),
            sum(r["seconds"][1] for r in results),
            sum(r["seconds"][2] for r in results),
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("data_dir")
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1, help="processes for the conversion")
    parser.add_argument("--check", type=int, default=200, help="joints per source to cross-check, about")
    parser.add_argument("--only", action="append", help="only the source with this name (repeatable)")
    parser.add_argument(
        "--smooth-within",
        type=float,
        default=0.0,
        metavar="DEGREES",
        help="count a joint as smooth if its kink is at most this in every master (default: exactly 0)",
    )
    args = parser.parse_args()

    paths = sorted(glob.glob(os.path.join(args.data_dir, "*.npz")))
    if args.only:
        paths = [p for p in paths if os.path.basename(p)[: -len(".npz")] in args.only]
    if not paths:
        raise SystemExit("no .npz files in %s" % args.data_dir)
    t0 = time.time()
    results = []
    pool = multiprocessing.Pool(args.jobs) if args.jobs > 1 else None
    try:
        for path in paths:
            t = time.time()
            results.append(process(path, pool, args))
            print("%s: %.1f s" % (results[-1]["name"], time.time() - t), file=sys.stderr, flush=True)
    finally:
        if pool:
            pool.close()
            pool.join()
    report(results, args, time.time() - t0, print)


if __name__ == "__main__":
    main()
