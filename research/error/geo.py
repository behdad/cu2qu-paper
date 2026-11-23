"""Vectorised building blocks: the spline family, its error polynomials, and the
pointwise ("normal-line") distance bound evaluated densely.  Used to separate
two questions: is the pointwise bound right and tight (this file, dense in u),
and what does a closed-form maximum over u lose (closed.py)."""
import numpy as np

def cross(u, v):
    return u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]

def dot(u, v):
    return (u * v).sum(-1)

def weights(mode, n):
    k = np.arange(n)
    if mode == 'ramp':
        return 1 - 2 * k / (n - 1) if n > 1 else np.zeros(1)
    if mode == 'zero':
        return np.zeros(n)
    if mode == 'alt':
        return np.where(k % 2 == 0, 1.0, -1.0)
    raise ValueError(mode)

def phi_coeffs(s):
    """Power-basis coefficients (n, 4) of phi_k(u) = c0 + c1 u + c2 u^2 + c3 u^3, and g (n+1)."""
    n = len(s)
    g = np.zeros(n + 1)
    g[1:n] = (s[:-1] + s[1:]) / 8
    h = (1 - s) / 2
    co = np.stack([g[:-1], -(2 * g[:-1] + h), 1 + g[:-1] + g[1:] + h, -np.ones(n)], -1)
    return co, g

def spline(c, s):
    """Off-curve (N, n, 2) and on-curve (N, n+1, 2) points of the spline with weights s."""
    n = len(s)
    c0, c1, c2, c3 = c[:, 0], c[:, 1], c[:, 2], c[:, 3]
    b1 = 3 * (c1 - c0); b2 = 3 * (c0 - 2 * c1 + c2); b3 = c3 - 3 * c2 + 3 * c1 - c0
    d = b3 / n ** 3
    off = np.empty((len(c), n, 2))
    for k in range(n):
        t0, t1 = k / n, (k + 1) / n
        B0 = c0 + t0 * (b1 + t0 * (b2 + t0 * b3)) if k else c0
        B1 = c0 + t1 * (b1 + t1 * (b2 + t1 * b3)) if k + 1 < n else c3
        D0 = b1 + t0 * (2 * b2 + 3 * t0 * b3)
        D1 = b1 + t1 * (2 * b2 + 3 * t1 * b3)
        p1 = B0 + D0 / (3 * n); p2 = B1 - D1 / (3 * n)
        qm = (3 * (p1 + p2) - B0 - B1) / 4
        off[:, k] = qm + s[k] * d / 4
    on = np.empty((len(c), n + 1, 2))
    on[:, 0] = c0; on[:, n] = c3
    on[:, 1:n] = (off[:, :-1] + off[:, 1:]) / 2
    return off, on, d

def pointwise_dense(c, s, U=129):
    """Dense evaluation of the pointwise bound for the spline with weights s (n >= 2).

    For piece k with legs a, z (V(u) = (1-u) a + u z, C = a x z), E = -d phi_k(u):
      first  = |E x V| / |V|                             error across the tangent
      eps    = (E.w) / (|V| (1 + sqrt(1 + y))),  y = ((z-a).w)(E.w) / |V|^2,  w = V/|V|
      exact  = |E.n - (C/|V|) eps^2|                     distance along the normal line
    Returns the maxima over all pieces and u of `first` and of `exact`, and a flag
    for cubics where y < -1 somewhere (no intersection with the normal line)."""
    n = len(s)
    off, on, d = spline(c, s)
    co, g = phi_coeffs(s)
    u = np.linspace(0, 1, U)[None, :]
    first_max = np.zeros(len(c)); exact_max = np.zeros(len(c)); bad = np.zeros(len(c), bool)
    for k in range(n):
        a = off[:, k] - on[:, k]; z = on[:, k + 1] - off[:, k]
        C = cross(a, z)
        phi = co[k, 0] + u * (co[k, 1] + u * (co[k, 2] + u * co[k, 3]))          # (1, U)
        V = (1 - u)[..., None] * a[:, None, :] + u[..., None] * z[:, None, :]   # (N, U, 2)
        Vn = np.sqrt(dot(V, V))
        with np.errstate(divide='ignore', invalid='ignore'):
            E = -d[:, None, :] * phi[..., None]
            Ew = dot(E, V) / Vn
            En = cross(V, E) / Vn                     # E . rot(w)
            y = dot((z - a)[:, None, :], V) / Vn * Ew / Vn ** 2
            eps = Ew / (Vn * (1 + np.sqrt(np.clip(1 + y, 0, None))))
            first = np.abs(En)
            exact = np.abs(En - C[:, None] / Vn * eps ** 2)
        degenerate = ~(Vn > 0)
        first = np.where(degenerate, np.sqrt(dot(E, E)), first)
        exact = np.where(degenerate, np.sqrt(dot(E, E)), exact)
        bad |= (np.nan_to_num(y, nan=0.0) < -1).any(1)
        first_max = np.maximum(first_max, first.max(1))
        exact_max = np.maximum(exact_max, exact.max(1))
    return first_max, exact_max, bad
