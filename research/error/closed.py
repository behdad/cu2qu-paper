"""Closed-form upper bounds on the geometric (Hausdorff) distance between a cubic
and its cu2qu spline.  Vectorised over cubics (numpy), O(1) arithmetic per piece,
no sampling / iteration / recursion.

  n1_closed(c)              n = 1, the tangent-intersection quadratic
  ramp_closed(c, n)         n >= 2, cu2qu's ramp
  spline_closed(c, n, s)    n >= 2, any weights s (uses the cubic formula, via
                            companion-matrix eigenvalues, to maximise a quartic)

Every function returns (bound, valid).  Where `valid` is False the conditions of
the proof do not hold and the caller must use the equal-parameter bound.
cu2qu_geo.py is the same thing for one cubic in plain Python.
"""
import math
import numpy as np
from geo import cross, dot, spline, phi_coeffs, weights

# ----------------------------------------------------------------------------
# small engines
# ----------------------------------------------------------------------------

def seg_min_norm(a, z):
    """min over u in [0,1] of |(1-u) a + u z|: distance from the origin to the segment."""
    w = z - a
    ww = dot(w, w)
    with np.errstate(divide='ignore', invalid='ignore'):
        u = np.clip(np.where(ww > 0, -dot(a, w) / ww, 0.0), 0, 1)
    v = a + u[..., None] * w
    return np.sqrt(dot(v, v))

def F_engine(al, be):
    """max over t in [0,1] of t^2 (1-t)^2 / (al (1-t) + be t), al, be > 0.
    The maximum is at the root in (0,1) of 3 (be-al) t^2 - (be - 5 al) t - 2 al = 0,
    written so that it stays accurate when be -> al."""
    disc = (be - 5 * al) ** 2 + 24 * al * (be - al)
    with np.errstate(divide='ignore', invalid='ignore'):
        t = 4 * al / ((5 * al - be) + np.sqrt(np.maximum(disc, 0)))
    t = np.clip(np.where(np.isfinite(t), t, 0.5), 0, 1)
    return t * t * (1 - t) ** 2 / (al * (1 - t) + be * t)

def H1_engine(al, be):
    """max over t in [0,1] of t (1-t)^2 / (al (1-t) + be t), al, be > 0:
    at t = 2 al / (3 al + sqrt(al (al + 8 be)))."""
    with np.errstate(divide='ignore', invalid='ignore'):
        t = 2 * al / (3 * al + np.sqrt(np.maximum(al * (al + 8 * be), 0)))
    t = np.clip(np.where(np.isfinite(t), t, 1 / 3), 0, 1)
    return t * (1 - t) ** 2 / (al * (1 - t) + be * t)

def quad_roots(qa, qb, qc):
    """The two roots of qa x^2 + qb x + qc (NaN where there is none), cancellation-free."""
    disc = qb * qb - 4 * qa * qc
    sq = np.sqrt(np.maximum(disc, 0))
    q = -0.5 * (qb + np.where(qb >= 0, 1.0, -1.0) * sq)
    with np.errstate(divide='ignore', invalid='ignore'):
        r1 = np.where(disc >= 0, q / qa, np.nan)
        r2 = np.where(disc >= 0, qc / q, np.nan)
    return r1, r2

def E_engine(A, p0, p1, argmax=False):
    """max over u in [0,1] of u^2 (A - u) max(p0 + p1 u, 0), A >= 1.
    Critical points: u = 1 and the roots of 4 p1 u^2 + 3 (p0 - A p1) u - 2 A p0 = 0."""
    def val(u):
        return u * u * (A - u) * np.maximum(p0 + p1 * u, 0)
    ubest = np.ones_like(p0)
    best = val(ubest)
    for u in quad_roots(4 * p1, 3 * (p0 - A * p1), -2 * A * p0):
        u = np.where(np.isfinite(u), np.clip(u, 0, 1), 1.0)
        v = val(u)
        better = v > best
        best = np.where(better, v, best); ubest = np.where(better, u, ubest)
    if argmax:
        return best, ubest
    return best

def G_engine(l0, l1, argmax=False):
    """max over t in [0,1] of t^2 (1-t)^2 max(l0 + l1 t, 0).
    Critical points: the roots of 5 l1 t^2 - (3 l1 - 4 l0) t - 2 l0 = 0."""
    def val(t):
        return t * t * (1 - t) ** 2 * np.maximum(l0 + l1 * t, 0)
    tbest = np.full_like(l0, 0.5)
    best = val(tbest)
    for t in quad_roots(5 * l1, -(3 * l1 - 4 * l0), -2 * l0):
        t = np.where(np.isfinite(t), np.clip(t, 0, 1), 0.5)
        v = val(t)
        better = v > best
        best = np.where(better, v, best); tbest = np.where(better, t, tbest)
    if argmax:
        return best, tbest
    return best

def mobius_majorant(n0, n1, d0, d1, ut):
    """A linear majorant p0 + p1 u on [0,1] of the Moebius function
    (n0 (1-u) + n1 u) / (d0 (1-u) + d1 u), d0, d1 > 0:
    its chord if it is convex, its tangent at ut if it is concave."""
    convex = (d1 - d0) * (n1 * d0 - n0 * d1) <= 0
    with np.errstate(divide='ignore', invalid='ignore'):
        c0 = n0 / d0; c1 = n1 / d1 - n0 / d0
        m = d0 * (1 - ut) + d1 * ut
        f = (n0 * (1 - ut) + n1 * ut) / m
        fp = (n1 * d0 - n0 * d1) / (m * m)
    return np.where(convex, c0, f - fp * ut), np.where(convex, c1, fp)

def mobius_envelopes(la, lz, al, be, ut):
    """Linear majorants on [0,1] of +mu and of -mu, mu = (la (1-u) + lz u) / (al (1-u) + be u),
    each with a flag saying whether that sign occurs on [0,1] at all."""
    out = []
    for sgn in (1.0, -1.0):
        p0, p1 = mobius_majorant(sgn * la, sgn * lz, al, be, ut)
        out.append((p0, p1, np.maximum(sgn * la, sgn * lz) > 0))
    return out

def quartic_max_abs(co, p0, p1):
    """max over u in [0,1] of |phi(u)| max(p0 + p1 u, 0), phi the cubic with power
    coefficients co (4,), p0, p1 arrays.  The critical points are the roots of a cubic;
    they are computed here as eigenvalues of its companion matrix."""
    N = len(p0)
    # f(u) = phi(u) (p0 + p1 u) = sum_j q_j u^j
    q = np.zeros((N, 5))
    for j in range(4):
        q[:, j] += co[j] * p0
        q[:, j + 1] += co[j] * p1
    dq = q[:, 1:] * np.arange(1, 5)            # cubic: dq0 + dq1 u + dq2 u^2 + dq3 u^3
    def val(u):
        ph = co[0] + u * (co[1] + u * (co[2] + u * co[3]))
        return np.abs(ph) * np.maximum(p0 + p1 * u, 0)
    best = np.maximum(val(np.zeros(N)), val(np.ones(N)))
    lead = dq[:, 3]
    scale = np.abs(dq).max(1)
    ok = np.abs(lead) > 1e-12 * scale
    comp = np.zeros((N, 3, 3))
    comp[:, 1, 0] = 1; comp[:, 2, 1] = 1
    with np.errstate(divide='ignore', invalid='ignore'):
        comp[:, 0, 2] = np.where(ok, -dq[:, 0] / lead, 0)
        comp[:, 1, 2] = np.where(ok, -dq[:, 1] / lead, 0)
        comp[:, 2, 2] = np.where(ok, -dq[:, 2] / lead, 0)
    ev = np.linalg.eigvals(comp)
    for j in range(3):
        u = np.clip(np.where(ok, ev[:, j].real, 0.5), 0, 1)
        best = np.maximum(best, val(u))
    # degenerate leading coefficient: the derivative is (at most) a quadratic
    for u in quad_roots(dq[:, 2], dq[:, 1], dq[:, 0]):
        u = np.where(np.isfinite(u) & ~ok, np.clip(u, 0, 1), 0.0)
        best = np.maximum(best, val(u))
    # zeros of the linear factor and a few fixed points cost nothing and guard the eigenvalue accuracy
    for u in (0.25, 0.5, 0.75):
        best = np.maximum(best, val(np.full(N, u)))
    return best

# ----------------------------------------------------------------------------
# n = 1
# ----------------------------------------------------------------------------

def n1_setup(c):
    """fontTools' tangent intersection q; legs a = q - c0, z = c3 - q; handle ratios."""
    c0, c1, c2, c3 = c[:, 0], c[:, 1], c[:, 2], c[:, 3]
    ab = c1 - c0; cd = c3 - c2
    pp = np.stack([-ab[:, 1], ab[:, 0]], -1)
    den = dot(pp, cd)
    den = np.where(np.abs(den) < 1e-15, 0.0, den)
    num = dot(pp, c0 - c2)
    num = np.where(np.abs(num) < 1e-15, 0.0, num)
    with np.errstate(divide='ignore', invalid='ignore'):
        h = np.where(den == 0, np.nan, num / den)
        q = c2 + cd * h[:, None]
        a = q - c0; z = c3 - q
        l1 = dot(c1 - c0, a) / dot(a, a)
        l2 = dot(c3 - c2, z) / dot(z, z)
    return q, a, z, l1, l2

def n1_closed(c, refined=True, detail=False):
    """Upper bound on the Hausdorff distance between a cubic and the quadratic whose
    control point q is the intersection of the cubic's end tangents.

      a = q - c0, z = c3 - q, c1 = c0 + l1 a, c2 = c3 - l2 z,  W = a - z,  K = |a x z|
      p = (3 l1 - 2) / 6,  r = (2 - 3 l2) / 6,  e(t) = (1-t) p + t r
      rho(t) = 3 t^2 (1-t)^2 g2(t),  g2(t) = 2 (p - r) + 3 e(t)^2     [B(t) = Q(tau(t)) + rho(t) W]

    Two bounds (the smaller valid one is returned):
      axis   = (3/16) max|g2| |W|                     needs only 0 <= l1, l2 <= 2
      normal = c 3 K max_t t^2 (1-t)^2 |g2(t)| / (al' (1-t) + be' t)
               where al' (1-t) + be' t is a linear minorant of |V(tau(t))|
    `refined` maximises the last expression with linear majorants (G_engine);
    otherwise it is bounded by max|g2| * F(al', be').
    """
    q, a, z, l1, l2 = n1_setup(c)
    N = len(c)
    with np.errstate(divide='ignore', invalid='ignore'):
        p = (3 * l1 - 2) / 6
        r = (2 - 3 * l2) / 6
        K = np.abs(cross(a, z))
        Wv = a - z
        Wn = np.sqrt(dot(Wv, Wv))
        A2 = 2 * (p - r)
        e2max = np.maximum(p * p, r * r)
        e2min = np.where(p * r <= 0, 0.0, np.minimum(p * p, r * r))
        hi = A2 + 3 * e2max                      # max of g2 on [0,1]
        lo = A2 + 3 * e2min                      # min of g2 on [0,1]
        M = np.maximum(np.abs(hi), np.abs(lo))
        axis = 3 / 16 * M * Wn
        axis_ok = np.isfinite(axis) & (l1 >= 0) & (l1 <= 2) & (l2 >= 0) & (l2 <= 2)

        s = a + z
        w = s / np.sqrt(dot(s, s))[:, None]
        al = dot(a, w); be = dot(z, w)
        m1 = al + (be - al) * l1 / 2
        m2 = al + (be - al) * (1 - l2 / 2)
        alp = np.minimum(al, 2 * m1 - m2)
        bep = np.minimum(be, 2 * m2 - m1)
        pos = (alp > 0) & (bep > 0)
        alp_ = np.where(pos, alp, 1.0); bep_ = np.where(pos, bep, 1.0)
        # points of the cubic outside the parabola (rho > 0): factor c
        xmax = 3 * np.maximum(hi, 0) * Wn ** 2 / (np.sqrt(alp_) + np.sqrt(bep_)) ** 4
        cfac = 2 / (1 + np.sqrt(np.clip(1 - xmax, 0, 1)))
        # the foot stays on the arc
        m0 = np.minimum(np.minimum(1.5 * l1, 1.5 - 0.75 * l2), 1.0)
        m1b = np.minimum(np.minimum(1.5 * l2, 1.5 - 0.75 * l1), 1.0)
        arc = (1.5 * cfac * M * Wn * H1_engine(alp_, bep_) <= m0) & (1.5 * cfac * M * Wn * H1_engine(bep_, alp_) <= m1b)
        simple = cfac * 3 * M * K * F_engine(alp_, bep_)
        normal = simple
        if refined:
            g0 = A2 + 3 * p * p; g1 = A2 + 3 * r * r            # g2 at t = 0, 1
            best = np.full(N, np.inf)
            tr = np.full(N, 0.5)
            for it in range(2):
                tot = np.zeros(N); tnew = tr.copy()
                # + part: g2 is convex, so below its chord
                l0, l1_ = mobius_majorant(g0, g1, alp_, bep_, tr)
                v, tt = G_engine(l0, l1_, argmax=True)
                v = np.where(hi > 0, v, 0.0)
                tnew = np.where(v > tot, tt, tnew); tot = np.maximum(tot, v)
                # - part: -g2 is concave, so below its tangent at tr
                er = (1 - tr) * p + tr * r
                gr = A2 + 3 * er * er; gpr = 6 * er * (r - p)
                n0 = -(gr - gpr * tr); n1 = -(gr + gpr * (1 - tr))
                l0, l1_ = mobius_majorant(n0, n1, alp_, bep_, tr)
                v, tt = G_engine(l0, l1_, argmax=True)
                v = np.where(lo < 0, v, 0.0)
                tnew = np.where(v > tot, tt, tnew); tot = np.maximum(tot, v)
                best = np.minimum(best, tot)
                tr = np.clip(tnew, 0.05, 0.95)
            normal = np.minimum(simple, cfac * 3 * K * best)
    normal_ok = (np.isfinite(normal) & (l1 > 0) & (l1 < 2) & (l2 > 0) & (l2 < 2) & pos & (xmax <= 0.75) & arc)
    bound = np.minimum(np.where(normal_ok, normal, np.inf), np.where(axis_ok, axis, np.inf))
    valid = normal_ok | axis_ok
    if detail:
        return bound, valid, dict(l1=l1, l2=l2, p=p, r=r, M=M, K=K, cfac=cfac, xmax=xmax, arc=arc, alp=alp, bep=bep,
                                  axis=axis, axis_ok=axis_ok, normal=normal, normal_ok=normal_ok, simple=simple, pos=pos)
    return bound, valid

# ----------------------------------------------------------------------------
# n >= 2
# ----------------------------------------------------------------------------

def piece_consts(s):
    """Per-piece constants of the weights s: phi coefficients, g, max|phi_k|, the jump J of
    phi'' at each joint, max|phi_k'|, and the parameter where |phi_k| is largest."""
    n = len(s)
    co, g = phi_coeffs(s)
    h = (1 - s) / 2
    Phi = np.zeros(n); D1 = np.zeros(n); upk = np.zeros(n)
    for k in range(n):
        c0_, c1_, c2_, c3_ = co[k]
        cand = [0.0, 1.0]
        disc = 4 * c2_ * c2_ - 12 * c3_ * c1_
        if disc >= 0:
            for sg in (1, -1):
                u = (-2 * c2_ + sg * math.sqrt(disc)) / (6 * c3_)
                if 0 < u < 1:
                    cand.append(u)
        vals = [abs(c0_ + u * (c1_ + u * (c2_ + u * c3_))) for u in cand]
        Phi[k] = max(vals); upk[k] = cand[int(np.argmax(vals))]
        cand = [abs(c1_), abs(c1_ + 2 * c2_ + 3 * c3_)]
        uv = -c2_ / (3 * c3_)
        if 0 < uv < 1:
            cand.append(abs(c1_ + 2 * c2_ * uv + 3 * c3_ * uv * uv))
        D1[k] = max(cand)
    J = np.zeros(n + 1)
    for k in range(1, n):
        J[k] = abs(6 + 2 * (g[k + 1] - g[k - 1] + h[k] - h[k - 1]))
    return co, g, Phi, J, D1, upk

def _piece_common(a, z, dn, Phi_k, D1_k, Jk, gk, extra=None):
    """T2, T3 and the validity of one piece.  `extra` = (A, beta) adds, for the first and
    last piece of the ramp, bounds that do not need the short leg."""
    C = np.abs(cross(a, z))
    Vmin = seg_min_norm(a, z)
    zma = z - a
    zn = np.sqrt(dot(zma, zma))
    with np.errstate(divide='ignore', invalid='ignore'):
        Ry = Phi_k / Vmin ** 2                       # sup |phi| / |V|^2
        R3 = Phi_k ** 2 / Vmin ** 3                  # sup phi^2 / |V|^3
        if extra is not None:
            A, beta = extra
            okb = beta > 0
            Ry = np.where(okb, np.minimum(Ry, A / beta ** 2), Ry)                 # phi_0 <= A u^2, |V| >= beta u
            R3 = np.where(okb, np.minimum(R3, 4 * A ** 3 / 27 / beta ** 3), R3)   # u (A-u)^2 <= 4 A^3 / 27
        Je = Jk if gk > 0 else 0.0                    # no foot crosses a joint where g = 0
        y0 = (zn + 0.5 * Je * dn) * dn * Ry
        root = 1 + np.sqrt(np.clip(1 - y0, 0, 1))
        T2 = C * dn * dn * R3 / root ** 2
        # |shift| <= |phi| epsc near the joints, and the shifted parameter stays in [0, 1]
        epsc = dn / (root * Vmin)
        okr = D1_k * epsc <= 1
        epsc = np.where(okr, epsc, np.inf)
        if extra is not None:
            eb = 2 * dn / (root * beta)              # |V| >= beta / 2 on the half next to the joint
            okrb = okb & (D1_k * eb <= 0.5)
            epsc = np.where(okrb, np.minimum(epsc, eb), epsc)
            okr = okr | okrb
        T3 = np.where(gk > 0, 0.5 * Jk * dn * (gk * np.where(okr, epsc, 0.0)) ** 2, 0.0)
        ok = (y0 <= 0.75) & okr & np.isfinite(T2) & np.isfinite(T3)
    return T2, T3, ok

def ramp_closed(c, n, refine=1, detail=False, engine='E'):
    """Upper bound on the Hausdorff distance between a cubic and its cu2qu spline of
    n >= 2 quadratics.  Per piece k, with legs a = Q_k - P_k, z = P_{k+1} - Q_k,
    V(u) = (1-u) a + u z, C = a x z, d = d3 / n^3:

      T1 = max_u |phi_k(u)| |d x V(u)| / |V(u)|       error across the tangent
      T2 = |C| |d|^2 sup(phi_k^2 / |V|^3) / root^2     curvature of the quadratic
      T3 = (J/2) |d| (max(|g_k|, |g_k+1|) |d| sup(1/|V|) / root)^2   change of quadratic at a joint
      root = 1 + sqrt(1 - y0),  y0 = (|z - a| + J |d| / 2) |d| sup(|phi_k| / |V|^2)

    T1 is bounded in closed form: for the first and the last piece by maximising
    u^2 (A - u) (p0 + p1 u) with p0 + p1 u a linear majorant of |d x V| / (V . w);
    for the pieces in between by Phi_k |d| max_u |sin(d, V(u))|.
    engine='prod' uses the latter (a product of two maxima) for every piece.
    """
    s = weights('ramp', n)
    co, g, Phi, J, D1, upk = piece_consts(s)
    A = 1 + g[1]                                    # phi_0 = u^2 (A - u)
    off, on, d = spline(c, s)
    dn = np.sqrt(dot(d, d))
    N = len(c)
    bound = np.zeros(N)
    valid = np.ones(N, bool)
    parts = np.zeros((N, n, 3))
    ustar = 2 * A / 3
    for k in range(n):
        a = off[:, k] - on[:, k]; z = on[:, k + 1] - off[:, k]
        Jk = max(J[k], J[k + 1]); gk = max(abs(g[k]), abs(g[k + 1]))
        if (k == 0 or k == n - 1) and engine == 'E':
            aa, zz = (a, z) if k == 0 else (z, a)           # the last piece is the first one mirrored
            la = cross(d, aa); lz = cross(d, zz)
            zero_a = (aa == 0).all(1)
            ut = np.full(N, ustar)
            T1 = None; beta0 = None
            for it in range(refine + 1):
                Vr = (1 - ut)[:, None] * aa + ut[:, None] * zz
                with np.errstate(divide='ignore', invalid='ignore'):
                    w = Vr / np.sqrt(dot(Vr, Vr))[:, None]
                al = dot(aa, w); be = dot(zz, w)
                if beta0 is None:
                    beta0 = np.where((al >= 0) & (be > 0), be, 0.0)
                good = (al > 0) & (be > 0)
                al_ = np.where(good, al, 1.0); be_ = np.where(good, be, 1.0)
                t1 = np.zeros(N); unew = ut.copy()
                for p0, p1, use in mobius_envelopes(la, lz, al_, be_, ut):
                    v, uu = E_engine(A, p0, p1, argmax=True)
                    v = np.where(use, v, 0.0)
                    unew = np.where(good & (v > t1), np.clip(uu, 0.02, 1.0), unew)   # next tangent point
                    t1 = np.maximum(t1, v)
                with np.errstate(divide='ignore', invalid='ignore'):
                    t1z = Phi[k] * np.abs(lz) / be          # a = 0: the tangent keeps the direction of z
                t1 = np.where(good, t1, np.where(zero_a & (be > 0), t1z, Phi[k] * dn))
                T1 = t1 if T1 is None else np.minimum(T1, t1)
                ut = unew
            T2, T3, ok = _piece_common(a, z, dn, Phi[k], D1[k], Jk, gk, extra=(A, beta0))
        else:
            da = dot(d, a); dz = dot(d, z)
            with np.errstate(divide='ignore', invalid='ignore'):
                sa = np.abs(cross(d, a)) / np.sqrt(dot(a, a))
                sz = np.abs(cross(d, z)) / np.sqrt(dot(z, z))
            smax = np.where(da * dz <= 0, dn, np.maximum(sa, sz))
            smax = np.where(np.isfinite(smax), smax, dn)
            T1 = Phi[k] * smax
            T2, T3, ok = _piece_common(a, z, dn, Phi[k], D1[k], Jk, gk)
        tot = T1 + T2 + T3
        parts[:, k, 0] = T1; parts[:, k, 1] = T2; parts[:, k, 2] = T3
        bound = np.maximum(bound, np.where(ok, tot, np.inf))
        valid &= ok
    if detail:
        return bound, valid, parts
    return bound, valid

def spline_closed(c, n, s, detail=False):
    """The same bound for arbitrary weights s (n >= 2): T1 of every piece is the exact
    maximum of |phi_k(u)| times a linear majorant of |d x V| / (V . w)."""
    co, g, Phi, J, D1, upk = piece_consts(s)
    off, on, d = spline(c, s)
    dn = np.sqrt(dot(d, d))
    N = len(c)
    bound = np.zeros(N)
    valid = np.ones(N, bool)
    parts = np.zeros((N, n, 3))
    for k in range(n):
        a = off[:, k] - on[:, k]; z = on[:, k + 1] - off[:, k]
        Jk = max(J[k], J[k + 1]); gk = max(abs(g[k]), abs(g[k + 1]))
        la = cross(d, a); lz = cross(d, z)
        ut = np.full(N, min(max(upk[k], 0.1), 0.9))
        Vr = (1 - ut)[:, None] * a + ut[:, None] * z
        with np.errstate(divide='ignore', invalid='ignore'):
            w = Vr / np.sqrt(dot(Vr, Vr))[:, None]
        al = dot(a, w); be = dot(z, w)
        good = (al > 0) & (be > 0)
        al_ = np.where(good, al, 1.0); be_ = np.where(good, be, 1.0)
        t1 = np.zeros(N)
        for p0, p1, use in mobius_envelopes(la, lz, al_, be_, ut):
            v = quartic_max_abs(co[k], p0, p1)
            t1 = np.maximum(t1, np.where(use, v, 0.0))
        T1 = np.where(good, t1, Phi[k] * dn)
        T2, T3, ok = _piece_common(a, z, dn, Phi[k], D1[k], Jk, gk)
        tot = T1 + T2 + T3
        parts[:, k, 0] = T1; parts[:, k, 1] = T2; parts[:, k, 2] = T3
        bound = np.maximum(bound, np.where(ok, tot, np.inf))
        valid &= ok
    if detail:
        return bound, valid, parts
    return bound, valid
