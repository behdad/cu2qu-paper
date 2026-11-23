#!/usr/bin/env python3
"""Geometric fit test for cu2qu splines: reference implementation, one cubic at a
time, plain Python (points are complex numbers, as in fontTools.cu2qu).

    n1_geo_bound(c0, c1, c2, c3)        upper bound on the Hausdorff distance between the cubic
                                        and its tangent-intersection quadratic, or None
    ramp_geo_bound(c0, c1, c2, c3, n)   the same for the cu2qu spline of n >= 2 quadratics, or None
    n1_ep_error / ramp_ep_error         the equal-parameter errors (what cu2qu tests today), exact
    quadratics_needed(c, tol)           least n whose bound (geometric, else equal-parameter) is <= tol

None means the conditions of the proof do not hold; the caller then uses the
equal-parameter error, which is always a valid bound.  Everything is O(1) per
quadratic: no sampling, no iteration to convergence, no recursion.  closed.py is
the numpy twin of this file; check_scalar.py compares the two.
"""
import math

INF = float('inf')

def dot(u, v):
    return u.real * v.real + u.imag * v.imag

def cross(u, v):
    return u.real * v.imag - u.imag * v.real

# ---------------------------------------------------------------------------
# engines: maxima of a few one-variable functions, each by one quadratic equation
# ---------------------------------------------------------------------------

def quad_roots(qa, qb, qc):
    """Real roots of qa x^2 + qb x + qc, without cancellation."""
    disc = qb * qb - 4 * qa * qc
    if disc < 0:
        return ()
    q = -0.5 * (qb + math.copysign(math.sqrt(disc), qb if qb != 0 else 1.0))
    out = []
    if qa != 0:
        out.append(q / qa)
    if q != 0:
        out.append(qc / q)
    return out

def F_engine(al, be):
    """max_t t^2 (1-t)^2 / (al (1-t) + be t), al, be > 0."""
    disc = (be - 5 * al) ** 2 + 24 * al * (be - al)
    t = 4 * al / ((5 * al - be) + math.sqrt(max(disc, 0.0)))
    t = min(max(t, 0.0), 1.0)
    return t * t * (1 - t) ** 2 / (al * (1 - t) + be * t)

def H1_engine(al, be):
    """max_t t (1-t)^2 / (al (1-t) + be t), al, be > 0."""
    t = 2 * al / (3 * al + math.sqrt(al * (al + 8 * be)))
    return t * (1 - t) ** 2 / (al * (1 - t) + be * t)

def E_engine(A, p0, p1):
    """(max, argmax) over u in [0,1] of u^2 (A - u) max(p0 + p1 u, 0)."""
    def val(u):
        return u * u * (A - u) * max(p0 + p1 * u, 0.0)
    best, ubest = val(1.0), 1.0
    for u in quad_roots(4 * p1, 3 * (p0 - A * p1), -2 * A * p0):
        if u == u and 0 < u < 1:
            v = val(u)
            if v > best:
                best, ubest = v, u
    return best, ubest

def G_engine(l0, l1):
    """(max, argmax) over t in [0,1] of t^2 (1-t)^2 max(l0 + l1 t, 0)."""
    def val(t):
        return t * t * (1 - t) ** 2 * max(l0 + l1 * t, 0.0)
    best, tbest = val(0.5), 0.5
    for t in quad_roots(5 * l1, -(3 * l1 - 4 * l0), -2 * l0):
        if t == t and 0 < t < 1:
            v = val(t)
            if v > best:
                best, tbest = v, t
    return best, tbest

def mobius_majorant(n0, n1, d0, d1, ut):
    """(p0, p1): p0 + p1 u >= (n0 (1-u) + n1 u) / (d0 (1-u) + d1 u) on [0,1]; d0, d1 > 0.
    The chord where the function is convex, its tangent at ut where it is concave."""
    if (d1 - d0) * (n1 * d0 - n0 * d1) <= 0:
        return n0 / d0, n1 / d1 - n0 / d0
    m = d0 * (1 - ut) + d1 * ut
    f = (n0 * (1 - ut) + n1 * ut) / m
    fp = (n1 * d0 - n0 * d1) / (m * m)
    return f - fp * ut, fp

def seg_min_norm(a, z):
    """min over u in [0,1] of |(1-u) a + u z|."""
    w = z - a
    ww = dot(w, w)
    u = min(max(-dot(a, w) / ww, 0.0), 1.0) if ww > 0 else 0.0
    return abs(a + u * w)

# ---------------------------------------------------------------------------
# n = 1
# ---------------------------------------------------------------------------

def tangent_intersection(c0, c1, c2, c3):
    """The quadratic's control point as fontTools computes it, or None."""
    ab = c1 - c0; cd = c3 - c2
    p = ab * 1j
    den = dot(p, cd)
    if abs(den) < 1e-15:
        return None
    num = dot(p, c0 - c2)
    if abs(num) < 1e-15:
        num = 0.0
    return c2 + cd * (num / den)

def n1_ep_error(c0, c1, c2, c3, q):
    """max_t |Q(t) - B(t)| exactly.  Q - B = 3 t (1-t) ((1-t) e1 + t e2); the interior
    critical points of its square are the roots of
        3 c t^3 + (5 b - 2 c) t^2 + (2 a - 3 b) t - a = 0,
    a = |e1|^2, b = e1.(e2 - e1), c = |e2 - e1|^2  (one cubic equation, solved in closed form)."""
    e1 = c0 + (q - c0) * (2 / 3) - c1
    e2 = c3 + (q - c3) * (2 / 3) - c2
    a = dot(e1, e1); b = dot(e1, e2 - e1); c = dot(e2 - e1, e2 - e1)
    def val(t):
        return 3 * t * (1 - t) * abs((1 - t) * e1 + t * e2)
    best = val(0.5)
    for t in cubic_roots(3 * c, 5 * b - 2 * c, 2 * a - 3 * b, -a):
        if 0 < t < 1:
            best = max(best, val(t))
    return best

def cubic_roots(a3, a2, a1, a0):
    """Real roots of a3 x^3 + a2 x^2 + a1 x + a0 (Cardano / Viete)."""
    scale = max(abs(a3), abs(a2), abs(a1), abs(a0))
    if scale == 0:
        return []
    if abs(a3) <= 1e-14 * scale:
        return list(quad_roots(a2, a1, a0)) if abs(a2) > 1e-14 * scale else ([-a0 / a1] if a1 != 0 else [])
    b = a2 / a3; c = a1 / a3; d = a0 / a3
    p = c - b * b / 3
    q = 2 * b ** 3 / 27 - b * c / 3 + d
    disc = q * q / 4 + p ** 3 / 27
    if disc > 0:
        s = math.sqrt(disc)
        u = math.copysign(abs(-q / 2 + s) ** (1 / 3), -q / 2 + s)
        v = math.copysign(abs(-q / 2 - s) ** (1 / 3), -q / 2 - s)
        return [u + v - b / 3]
    if p == 0:
        return [-b / 3]
    m = 2 * math.sqrt(-p / 3)
    th = math.acos(max(-1.0, min(1.0, 3 * q / (p * m)))) / 3
    return [m * math.cos(th - 2 * math.pi * k / 3) - b / 3 for k in range(3)]

def n1_geo_bound(c0, c1, c2, c3, q=None, refined=True):
    """Upper bound on the Hausdorff distance between the cubic and the quadratic (c0, q, c3),
    q the intersection of the end tangents.  None if the conditions do not hold."""
    if q is None:
        q = tangent_intersection(c0, c1, c2, c3)
        if q is None:
            return None
    a = q - c0; z = c3 - q
    aa = dot(a, a); zz = dot(z, z)
    if not (aa > 0 and zz > 0):
        return None
    l1 = dot(c1 - c0, a) / aa               # c1 = c0 + l1 a
    l2 = dot(c3 - c2, z) / zz               # c2 = c3 - l2 z
    p = (3 * l1 - 2) / 6
    r = (2 - 3 * l2) / 6
    K = abs(cross(a, z))
    Wn = abs(a - z)
    A2 = 2 * (p - r)
    hi = A2 + 3 * max(p * p, r * r)                             # max of g2(t) = A2 + 3 ((1-t) p + t r)^2
    lo = A2 + 3 * (0.0 if p * r <= 0 else min(p * p, r * r))    # min
    M = max(abs(hi), abs(lo))
    best = INF
    # (1) along the axis of the parabola: B(t) = Q(tau(t)) + rho(t) (a - z), |rho| <= 3 M / 16
    if 0 <= l1 <= 2 and 0 <= l2 <= 2:
        best = 3 / 16 * M * Wn
    # (2) across the parabola: |rho| |a x z| / |V(tau)|, V(tau) = (1-tau) a + tau z
    s = a + z
    sn = abs(s)
    if sn > 0 and 0 < l1 < 2 and 0 < l2 < 2:
        w = s / sn
        al = dot(a, w); be = dot(z, w)
        m1 = al + (be - al) * l1 / 2
        m2 = al + (be - al) * (1 - l2 / 2)
        alp = min(al, 2 * m1 - m2)          # al' (1-t) + be' t <= V(tau(t)) . w <= |V(tau(t))|
        bep = min(be, 2 * m2 - m1)
        if alp > 0 and bep > 0:
            x = 3 * max(hi, 0.0) * Wn * Wn / (math.sqrt(alp) + math.sqrt(bep)) ** 4
            if x <= 0.75:
                cf = 2 / (1 + math.sqrt(1 - x))
                m0 = min(1.5 * l1, 1.5 - 0.75 * l2, 1.0)
                m1b = min(1.5 * l2, 1.5 - 0.75 * l1, 1.0)
                if (1.5 * cf * M * Wn * H1_engine(alp, bep) <= m0 and
                        1.5 * cf * M * Wn * H1_engine(bep, alp) <= m1b):
                    normal = cf * 3 * M * K * F_engine(alp, bep)
                    if refined:
                        g0 = A2 + 3 * p * p; g1 = A2 + 3 * r * r
                        tr = 0.5
                        fine = INF
                        for _ in range(2):
                            tot, tnew = 0.0, tr
                            if hi > 0:          # g2 is convex: below its chord
                                l0, l1_ = mobius_majorant(g0, g1, alp, bep, tr)
                                v, tt = G_engine(l0, l1_)
                                if v > tot:
                                    tot, tnew = v, tt
                            if lo < 0:          # -g2 is concave: below its tangent at tr
                                er = (1 - tr) * p + tr * r
                                gr = A2 + 3 * er * er; gpr = 6 * er * (r - p)
                                l0, l1_ = mobius_majorant(-(gr - gpr * tr), -(gr + gpr * (1 - tr)), alp, bep, tr)
                                v, tt = G_engine(l0, l1_)
                                if v > tot:
                                    tot, tnew = v, tt
                            fine = min(fine, tot)
                            tr = min(max(tnew, 0.05), 0.95)
                        normal = min(normal, cf * 3 * K * fine)
                    best = min(best, normal)
    return best if best < INF else None

# ---------------------------------------------------------------------------
# n >= 2: cu2qu's ramp
# ---------------------------------------------------------------------------

def ramp_consts(n):
    """Constants of the ramp that depend on n only: A, g_k, Phi_k, J_k, D1_k."""
    s = [1 - 2 * k / (n - 1) for k in range(n)]
    g = [0.0] + [(s[k - 1] + s[k]) / 8 for k in range(1, n)] + [0.0]
    h = [(1 - x) / 2 for x in s]
    A = 1 + g[1]
    cc = (n - 3) / (4.0 * (n - 1))
    Phi = [0.0] * n; D1 = [0.0] * n
    for k in range(n):
        # phi_k(u) = c0 + c1 u + c2 u^2 - u^3
        c0_, c1_, c2_ = g[k], -(2 * g[k] + h[k]), 1 + g[k] + g[k + 1] + h[k]
        cand = [0.0, 1.0] + [u for u in quad_roots(-3.0, 2 * c2_, c1_) if 0 < u < 1]
        Phi[k] = max(abs(c0_ + u * (c1_ + u * (c2_ - u))) for u in cand)
        cand = [0.0, 1.0] + ([c2_ / 3] if 0 < c2_ / 3 < 1 else [])
        D1[k] = max(abs(c1_ + 2 * c2_ * u - 3 * u * u) for u in cand)
    J = [0.0] * (n + 1)
    for k in range(1, n):
        J[k] = abs(6 + 2 * (g[k + 1] - g[k - 1] + h[k] - h[k - 1]))
    return s, g, A, Phi, J, D1

_consts = {}

def ramp_ep_error(d3n, n):
    """Equal-parameter error of the ramp spline: |d3| Phi(n) / n^3."""
    if n not in _consts:
        _consts[n] = ramp_consts(n)
    return d3n * max(_consts[n][3]) / n ** 3

def ramp_spline(c0, c1, c2, c3, n):
    """Off-curve points of the cu2qu spline of n >= 2 quadratics, and d = d3 / n^3."""
    b1 = 3 * (c1 - c0); b2 = 3 * (c0 - 2 * c1 + c2); b3 = c3 - 3 * c2 + 3 * c1 - c0
    d = b3 / n ** 3
    off = []
    for k in range(n):
        t0, t1 = k / n, (k + 1) / n
        B0 = c0 + t0 * (b1 + t0 * (b2 + t0 * b3)) if k else c0
        B1 = c0 + t1 * (b1 + t1 * (b2 + t1 * b3)) if k + 1 < n else c3
        D0 = b1 + t0 * (2 * b2 + 3 * t0 * b3)
        D1 = b1 + t1 * (2 * b2 + 3 * t1 * b3)
        p1 = B0 + D0 / (3 * n); p2 = B1 - D1 / (3 * n)
        qm = (3 * (p1 + p2) - B0 - B1) / 4
        off.append(qm + (1 - 2 * k / (n - 1)) * d / 4)
    return off, d

def _piece_terms(a, z, dn, Phi_k, D1_k, Jk, gk, A=None, beta=0.0):
    """(T2, T3) of one piece, or None if the piece fails the conditions."""
    C = abs(cross(a, z))
    Vmin = seg_min_norm(a, z)
    zn = abs(z - a)
    Ry = Phi_k / Vmin ** 2 if Vmin > 0 else INF            # sup |phi| / |V|^2
    R3 = Phi_k ** 2 / Vmin ** 3 if Vmin > 0 else INF       # sup phi^2 / |V|^3
    if beta > 0:                                           # first / last piece: |V(u)| >= beta u
        Ry = min(Ry, A / beta ** 2)
        R3 = min(R3, 4 * A ** 3 / 27 / beta ** 3)
    if dn == 0:
        return 0.0, 0.0
    if Ry == INF:
        return None
    Je = Jk if gk > 0 else 0.0
    y0 = (zn + 0.5 * Je * dn) * dn * Ry
    if not y0 <= 0.75:
        return None
    root = 1 + math.sqrt(1 - y0)
    T2 = C * dn * dn * R3 / root ** 2
    epsc = INF
    if Vmin > 0 and D1_k * dn / (root * Vmin) <= 1:
        epsc = dn / (root * Vmin)
    if beta > 0 and D1_k * 2 * dn / (root * beta) <= 0.5:
        epsc = min(epsc, 2 * dn / (root * beta))
    if epsc == INF:
        return None
    T3 = 0.5 * Jk * dn * (gk * epsc) ** 2 if gk > 0 else 0.0
    return T2, T3

def ramp_geo_bound(c0, c1, c2, c3, n, off=None, refine=1):
    """Upper bound on the Hausdorff distance between the cubic and its cu2qu spline of n >= 2
    quadratics.  None if the conditions do not hold for some piece."""
    if n not in _consts:
        _consts[n] = ramp_consts(n)
    s, g, A, Phi, J, D1 = _consts[n]
    if off is None:
        off, d = ramp_spline(c0, c1, c2, c3, n)
    else:
        d = (c3 - 3 * c2 + 3 * c1 - c0) / n ** 3
    dn = abs(d)
    on = [c0] + [(off[k - 1] + off[k]) / 2 for k in range(1, n)] + [c3]
    bound = 0.0
    for k in range(n):
        a = off[k] - on[k]; z = on[k + 1] - off[k]
        Jk = max(J[k], J[k + 1]); gk = max(abs(g[k]), abs(g[k + 1]))
        beta0 = 0.0
        if k == 0 or k == n - 1:
            aa, zz = (a, z) if k == 0 else (z, a)       # the last piece is the first one mirrored
            la = cross(d, aa); lz = cross(d, zz)
            ut = 2 * A / 3
            T1 = INF
            for it in range(refine + 1):
                Vr = (1 - ut) * aa + ut * zz
                vn = abs(Vr)
                if not vn > 0:
                    T1 = min(T1, Phi[k] * dn); break
                w = Vr / vn
                al = dot(aa, w); be = dot(zz, w)
                if it == 0 and al >= 0 and be > 0:
                    beta0 = be
                if al > 0 and be > 0:
                    t1, unew = 0.0, ut
                    for sgn in (1.0, -1.0):
                        if max(sgn * la, sgn * lz) > 0:
                            p0, p1 = mobius_majorant(sgn * la, sgn * lz, al, be, ut)
                            v, uu = E_engine(A, p0, p1)
                            if v > t1:
                                t1, unew = v, min(max(uu, 0.02), 1.0)
                    T1 = min(T1, t1)
                    ut = unew
                elif aa == 0 and be > 0:
                    T1 = min(T1, Phi[k] * abs(lz) / be)     # a = 0: the tangent keeps the direction of z
                else:
                    T1 = min(T1, Phi[k] * dn)
        else:
            an = abs(a); zn_ = abs(z)
            if dot(d, a) * dot(d, z) <= 0 or an == 0 or zn_ == 0:
                smax = dn
            else:
                smax = max(abs(cross(d, a)) / an, abs(cross(d, z)) / zn_)
            T1 = Phi[k] * smax
        terms = _piece_terms(a, z, dn, Phi[k], D1[k], Jk, gk, A, beta0)
        if terms is None:
            return None
        bound = max(bound, T1 + terms[0] + terms[1])
    return bound

# ---------------------------------------------------------------------------
# the rule for n
# ---------------------------------------------------------------------------

def bounds(c, n):
    """(equal-parameter error, geometric bound or None) for n quadratics; (inf, None) if n = 1 has no quadratic."""
    c0, c1, c2, c3 = c
    if n == 1:
        q = tangent_intersection(c0, c1, c2, c3)
        if q is None:
            return INF, None
        return n1_ep_error(c0, c1, c2, c3, q), n1_geo_bound(c0, c1, c2, c3, q)
    return ramp_ep_error(abs(c3 - 3 * c2 + 3 * c1 - c0), n), ramp_geo_bound(c0, c1, c2, c3, n)

def quadratics_needed(c, tol, max_n=100):
    """Least n whose spline is certified within tol of the cubic (as point sets)."""
    for n in range(1, max_n + 1):
        ep, geo = bounds(c, n)
        if ep <= tol or (geo is not None and geo <= tol):
            return n
    return None

if __name__ == '__main__':
    import sys
    c = [complex(232, 367), complex(255, 440), complex(276, 512), complex(297, 588)]
    for n in range(1, 5):
        print(n, bounds(c, n))
    print(quadratics_needed(c, 1.0))
