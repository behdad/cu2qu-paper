// REIMPLEMENTATIONS.  Nothing in this file is the original authors' code; each
// function is written from the source or the paper named in its comment, and
// is labelled "reimpl_*" in every table.

use crate::Seg;
use kurbo::{CubicBez, ParamCurve, Point, Vec2};

fn d3(c: &CubicBez) -> Vec2 {
    c.p3.to_vec2() - 3.0 * c.p2.to_vec2() + 3.0 * c.p1.to_vec2() - c.p0.to_vec2()
}

fn piece(c: &CubicBez, k: usize, n: usize) -> CubicBez {
    if n == 1 {
        return *c;
    }
    let mut s = c.subsegment(k as f64 / n as f64..(k + 1) as f64 / n as f64);
    if k == 0 {
        s.p0 = c.p0;
    }
    if k + 1 == n {
        s.p3 = c.p3;
    }
    s
}

/// The quadratic with the "mid-point" control point (Colomitchi 2006, config. 3;
/// kurbo ToQuads; lyon to_quadratic; Skia SkDCubic::toQuad).
fn qm(c: &CubicBez) -> Seg {
    let ctrl = (3.0 * (c.p1.to_vec2() + c.p2.to_vec2()) - c.p0.to_vec2() - c.p3.to_vec2()) / 4.0;
    Seg::quad(c.p0, ctrl.to_point(), c.p3)
}

fn start_tangent_point(c: &CubicBez) -> Point {
    c.p0 + (c.p1 - c.p0) * 1.5
}

fn end_tangent_point(c: &CubicBez) -> Point {
    c.p3 + (c.p2 - c.p3) * 1.5
}

/// Turn a list of off-curve points into explicit quads with implied on-curve points.
fn spline_from_offs(p0: Point, offs: &[Point], p3: Point) -> Vec<Seg> {
    let n = offs.len();
    let mut out = Vec::with_capacity(n);
    let mut from = p0;
    for k in 0..n {
        let to = if k + 1 < n { offs[k].midpoint(offs[k + 1]) } else { p3 };
        out.push(Seg::quad(from, offs[k], to));
        from = to;
    }
    out
}

/// General family: off-curve point of piece k = q_M,k + s_k d / 4, i.e. the
/// point at parameter (1 - s_k)/2 between the piece's start-tangent point (s = +1)
/// and its end-tangent point (s = -1).
pub fn family_spline(c: &CubicBez, s: &[f64]) -> Vec<Seg> {
    let n = s.len();
    let offs: Vec<Point> = (0..n)
        .map(|k| {
            let p = piece(c, k, n);
            start_tangent_point(&p).lerp(end_tangent_point(&p), (1.0 - s[k]) / 2.0)
        })
        .collect();
    spline_from_offs(c.p0, &offs, c.p3)
}

/// fontTools cu2qu cubic_approx_spline construction for n >= 2 (no error test).
pub fn ramp_spline(c: &CubicBez, n: usize) -> Vec<Seg> {
    let offs: Vec<Point> = (0..n)
        .map(|k| {
            let p = piece(c, k, n);
            start_tangent_point(&p).lerp(end_tangent_point(&p), k as f64 / (n - 1) as f64)
        })
        .collect();
    spline_from_offs(c.p0, &offs, c.p3)
}

/// fontTools cu2qu cubic_farthest_fit_inside.
fn fit_inside(p0: Vec2, p1: Vec2, p2: Vec2, p3: Vec2, tol: f64, depth: u32) -> bool {
    if p2.hypot() <= tol && p1.hypot() <= tol {
        return true;
    }
    if depth > 40 {
        return false;
    }
    let mid = (p0 + 3.0 * (p1 + p2) + p3) * 0.125;
    if mid.hypot() > tol {
        return false;
    }
    let deriv3 = (p3 + p2 - p1 - p0) * 0.125;
    fit_inside(p0, (p0 + p1) * 0.5, mid - deriv3, mid, tol, depth + 1)
        && fit_inside(mid, mid + deriv3, (p2 + p3) * 0.5, p3, tol, depth + 1)
}

fn cu2qu_dot(a: Vec2, b: Vec2) -> f64 {
    let r = a.x * b.x + a.y * b.y;
    if r.abs() < 1e-15 {
        0.0
    } else {
        r
    }
}

/// fontTools cu2qu cubic_approx_quadratic: the tangent-intersection quadratic, if it fits.
fn cu2qu_single(c: &CubicBez, tol: f64) -> Option<Point> {
    let ab = c.p1 - c.p0;
    let cd = c.p3 - c.p2;
    let p = Vec2::new(-ab.y, ab.x);
    let den = cu2qu_dot(p, cd);
    let q1 = if den == 0.0 {
        if c.p1 == c.p2 && (c.p0 == c.p1 || c.p2 == c.p3) {
            c.p1
        } else {
            return None;
        }
    } else {
        let h = cu2qu_dot(p, c.p0 - c.p2) / den;
        c.p2 + cd * h
    };
    let c1 = c.p0 + (q1 - c.p0) * (2.0 / 3.0);
    let c2 = c.p3 + (q1 - c.p3) * (2.0 / 3.0);
    if fit_inside(Vec2::ZERO, c1 - c.p1, c2 - c.p2, Vec2::ZERO, tol, 0) {
        Some(q1)
    } else {
        None
    }
}

/// The tangent-intersection quadratic as cu2qu computes it, without the error test.
pub fn tangent_single(c: &CubicBez) -> Option<Vec<Seg>> {
    let ab = c.p1 - c.p0;
    let cd = c.p3 - c.p2;
    let p = Vec2::new(-ab.y, ab.x);
    let den = cu2qu_dot(p, cd);
    let q1 = if den == 0.0 {
        if c.p1 == c.p2 && (c.p0 == c.p1 || c.p2 == c.p3) {
            c.p1
        } else {
            return None;
        }
    } else {
        let h = cu2qu_dot(p, c.p0 - c.p2) / den;
        c.p2 + cd * h
    };
    Some(vec![Seg::quad(c.p0, q1, c.p3)])
}

/// n equal pieces, each piece's control point at the intersection of the tangents
/// at its two ends (Groleau 2002 "tangent approach", FontForge's old method,
/// cubic2quad's processSegment), on-curve points on the cubic.  None if some
/// piece's tangents do not intersect.
pub fn tangent_uniform(c: &CubicBez, n: usize) -> Option<Vec<Seg>> {
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let p = piece(c, k, n);
        let mut ab = p.p1 - p.p0;
        if ab.hypot2() == 0.0 {
            ab = p.p2 - p.p0;
        }
        let mut cd = p.p3 - p.p2;
        if cd.hypot2() == 0.0 {
            cd = p.p3 - p.p1;
        }
        let den = ab.cross(cd);
        if den.abs() <= 1e-9 * ab.hypot() * cd.hypot() {
            return None;
        }
        let h = ab.cross(p.p0 - p.p3) / den;
        // intersection, written from the end point along the end tangent
        let q = p.p3 + cd * h;
        out.push(Seg::quad(p.p0, q, p.p3));
    }
    Some(out)
}

/// The same off-curve points as tangent_uniform, but with the on-curve points
/// implied (midpoints), TrueType-style.  For n = 2 this is what Roboto's 2014
/// convertCurves.py did ("a quadratic spline with exactly two off curve points").
pub fn tangent_implied(c: &CubicBez, n: usize) -> Option<Vec<Seg>> {
    let t = tangent_uniform(c, n)?;
    let offs: Vec<Point> = t.iter().map(|s| s.p1).collect();
    Some(spline_from_offs(c.p0, &offs, c.p3))
}

/// Skia's choice of control points for n equal pieces of an inflection-free cubic:
/// first piece on its start tangent, last piece on its end tangent, mid-point in
/// between (one piece: mid-point); on-curve points on the cubic.
pub fn skia_pattern(c: &CubicBez, n: usize) -> Vec<Seg> {
    (0..n)
        .map(|k| {
            let p = piece(c, k, n);
            let ctrl = if n == 1 {
                qm(&p).p1
            } else if k == 0 {
                start_tangent_point(&p)
            } else if k + 1 == n {
                end_tangent_point(&p)
            } else {
                qm(&p).p1
            };
            Seg::quad(p.p0, ctrl, p.p3)
        })
        .collect()
}

/// HarfBuzz hb_gpu_cubic_spline_error: equal-parameter error of the cu2qu spline
/// of n >= 2 quadratics, as a multiple of |d3| / n^3.
pub fn phi(n: usize) -> f64 {
    if n == 2 {
        return 4.0 / 27.0;
    }
    let nf = n as f64;
    let g = (nf - 2.0) / (4.0 * (nf - 1.0));
    let first = 4.0 / 27.0 * (1.0 + g).powi(3);
    let w2 = (0.25 - 0.5 / (nf - 1.0)) / 3.0;
    let second = (nf - 3.0) / (4.0 * (nf - 1.0)) + 2.0 * w2 * w2.sqrt();
    first.max(second)
}

/// cu2qu's spline with n computed directly from Phi(n) (what HarfBuzz's
/// hb-gpu-cu2qu.hh does), n = 1 decided by cu2qu's own single-quadratic test.
pub fn hb_direct(c: &CubicBez, tol: f64) -> Vec<Seg> {
    if let Some(q1) = cu2qu_single(c, tol) {
        return vec![Seg::quad(c.p0, q1, c.p3)];
    }
    let d = d3(c).hypot();
    let mut n = 2;
    while n <= 100 && !(d * phi(n) <= tol * (n * n * n) as f64) {
        n += 1;
    }
    if n > 100 {
        return vec![];
    }
    ramp_spline(c, n)
}

/// cu2qu's trial loop with cu2qu's own per-piece test (cubic_approx_spline).
pub fn cu2qu_trial(c: &CubicBez, tol: f64) -> Vec<Seg> {
    if let Some(q1) = cu2qu_single(c, tol) {
        return vec![Seg::quad(c.p0, q1, c.p3)];
    }
    for n in 2..=100 {
        let sp = ramp_spline(c, n);
        let mut ok = true;
        let mut d0 = Vec2::ZERO;
        for k in 0..n {
            let p = piece(c, k, n);
            let q = &sp[k];
            let d1 = q.p2 - p.p3;
            if d1.hypot() > tol
                || !fit_inside(
                    d0,
                    q.p0.lerp(q.p1, 2.0 / 3.0) - p.p1,
                    q.p2.lerp(q.p1, 2.0 / 3.0) - p.p2,
                    d1,
                    tol,
                    0,
                )
            {
                ok = false;
                break;
            }
            d0 = d1;
        }
        if ok {
            return sp;
        }
    }
    vec![]
}

const S3_36: f64 = 0.048_112_522_432_468_81; // sqrt(3) / 36

/// Colomitchi 2006, "Quadratic approximation of a cubic Bezier by half-splitting",
/// with the mid-point approximation and its defect sqrt(3)/36 |d3|.
pub fn colomitchi_half(c: &CubicBez, tol: f64) -> Vec<Seg> {
    fn rec(c: &CubicBez, tol: f64, depth: u32, out: &mut Vec<Seg>) {
        if S3_36 * d3(c).hypot() <= tol || depth >= 16 {
            out.push(qm(c));
        } else {
            let (l, r) = c.subdivide();
            rec(&l, tol, depth + 1, out);
            rec(&r, tol, depth + 1, out);
        }
    }
    let mut out = Vec::new();
    rec(c, tol, 0, &mut out);
    out
}

/// Colomitchi 2006, "simple adaptive division using mid-point approximation":
/// t_max from sqrt(3)/18 D01 t_max^3 = prec, cut there, repeat on the rest.
pub fn colomitchi_adaptive(c: &CubicBez, tol: f64) -> Vec<Seg> {
    let mut out = Vec::new();
    let mut cur = *c;
    for _ in 0..1000 {
        let d = S3_36 * d3(&cur).hypot();
        let tmax = if d > 0.0 { (tol / d).cbrt() } else { f64::INFINITY };
        if tmax >= 1.0 {
            out.push(qm(&cur));
            break;
        }
        let mut first = cur.subsegment(0.0..tmax);
        first.p0 = cur.p0;
        out.push(qm(&first));
        let mut rest = cur.subsegment(tmax..1.0);
        rest.p0 = first.p3;
        rest.p3 = cur.p3;
        cur = rest;
    }
    out
}

/// Colomitchi 2006, "adaptive symmetric division": t_max and 1 - t_max at once.
pub fn colomitchi_symmetric(c: &CubicBez, tol: f64) -> Vec<Seg> {
    let mut left = Vec::new();
    let mut right = Vec::new();
    let mut cur = *c;
    for _ in 0..1000 {
        let d = S3_36 * d3(&cur).hypot();
        let tmax = if d > 0.0 { (tol / d).cbrt() } else { f64::INFINITY };
        if tmax >= 1.0 {
            left.push(qm(&cur));
            break;
        }
        if tmax >= 0.5 {
            let (l, r) = cur.subdivide();
            left.push(qm(&l));
            left.push(qm(&r));
            break;
        }
        let mut a = cur.subsegment(0.0..tmax);
        a.p0 = cur.p0;
        let mut b = cur.subsegment(1.0 - tmax..1.0);
        b.p3 = cur.p3;
        let mut mid = cur.subsegment(tmax..1.0 - tmax);
        mid.p0 = a.p3;
        mid.p3 = b.p0;
        left.push(qm(&a));
        right.push(qm(&b));
        cur = mid;
    }
    right.reverse();
    left.extend(right);
    left
}

/// Truong, Yuksel, Seiler 2020 (eqs. 13-15, gamma = 1/2): two quadratics per cubic
/// piece.  Here the cubic is first cut into m equal pieces ("pre-subdividing",
/// their section 4), with m the least for which eq. 29's error |a3| / 54 / m^3 is
/// within tolerance.  The paper states the error and suggests subdividing by it;
/// the closed form for m is this file's, not the paper's.  n = 1 is decided by
/// cu2qu's single-quadratic test so that counts can be compared with cu2qu.
pub fn truong_presub(c: &CubicBez, tol: f64) -> Vec<Seg> {
    if let Some(q1) = cu2qu_single(c, tol) {
        return vec![Seg::quad(c.p0, q1, c.p3)];
    }
    let d = d3(c).hypot();
    let m = ((d / (54.0 * tol)).cbrt().ceil() as usize).max(1);
    truong_m(c, m)
}

pub fn truong_m(c: &CubicBez, m: usize) -> Vec<Seg> {
    let mut offs = Vec::with_capacity(2 * m);
    for k in 0..m {
        let p = piece(c, k, m);
        offs.push(p.p0 + (p.p1 - p.p0) * 0.75);
        offs.push(p.p3 + (p.p2 - p.p3) * 0.75);
    }
    spline_from_offs(c.p0, &offs, c.p3)
}

pub fn truong2_fixed(c: &CubicBez, _tol: f64) -> Vec<Seg> {
    truong_m(c, 1)
}

/// Groleau 2002, drawCubicBezier2 ("fixed MidPoint approach", after Helen Triolo):
/// always four quadratics.  Transcribed from bezier_lib.as.
pub fn triolo4_fixed(c: &CubicBez, _tol: f64) -> Vec<Seg> {
    let pa = c.p0.lerp(c.p1, 0.75);
    let pb = c.p3.lerp(c.p2, 0.75);
    let d = (c.p3 - c.p0) / 16.0;
    let pc1 = c.p0.lerp(c.p1, 0.375);
    let pc2 = pa.lerp(pb, 0.375) - d;
    let pc3 = pb.lerp(pa, 0.375) + d;
    let pc4 = c.p3.lerp(c.p2, 0.375);
    let a1 = pc1.midpoint(pc2);
    let a2 = pa.midpoint(pb);
    let a3 = pc3.midpoint(pc4);
    vec![
        Seg::quad(c.p0, pc1, a1),
        Seg::quad(a1, pc2, a2),
        Seg::quad(a2, pc3, a3),
        Seg::quad(a3, pc4, c.p3),
    ]
}

/// Groleau 2002, "BSpline approach": two quadratics with control points P1 and P2
/// meeting at the middle of [P1, P2].
pub fn groleau_bspline2_fixed(c: &CubicBez, _tol: f64) -> Vec<Seg> {
    spline_from_offs(c.p0, &[c.p1, c.p2], c.p3)
}

/// Alfke 1994 / Penner's "generic MidPoint approach" in Groleau 2002: the quadratic
/// whose control point is the intersection of the end tangents; the metric is the
/// distance between the midpoints of the quadratic and of the cubic; halve and
/// recurse while it exceeds the tolerance.  (Parallel tangents: neither source
/// says; the mid-point control point is used here.)
pub fn alfke_penner(c: &CubicBez, tol: f64) -> Vec<Seg> {
    fn rec(c: &CubicBez, tol: f64, depth: u32, out: &mut Vec<Seg>) {
        let ab = c.p1 - c.p0;
        let cd = c.p3 - c.p2;
        let den = ab.cross(cd);
        let e = if den != 0.0 {
            let h = ab.cross(c.p0 - c.p2) / den;
            c.p2 + cd * h
        } else {
            qm(c).p1
        };
        let qmid = (c.p0.to_vec2() + 2.0 * e.to_vec2() + c.p3.to_vec2()) / 4.0;
        let cmid = (c.p0.to_vec2() + 3.0 * (c.p1.to_vec2() + c.p2.to_vec2()) + c.p3.to_vec2()) / 8.0;
        if (qmid - cmid).hypot() <= tol || depth >= 10 {
            out.push(Seg::quad(c.p0, e, c.p3));
        } else {
            let (l, r) = c.subdivide();
            rec(&l, tol, depth + 1, out);
            rec(&r, tol, depth + 1, out);
        }
    }
    let mut out = Vec::new();
    rec(c, tol, 0, &mut out);
    out
}

// ------------------------------------------------------------------ constants

fn eq_const(c: &CubicBez, segs: &[Seg]) -> f64 {
    // max over t of |C(t) - S(t)| at equal parameters, over |d3| / n^3
    let n = segs.len();
    let mut best: f64 = 0.0;
    for k in 0..n {
        let q = kurbo::QuadBez::new(segs[k].p0, segs[k].p1, segs[k].p2);
        for j in 0..=20000 {
            let u = j as f64 / 20000.0;
            let d = (c.eval((k as f64 + u) / n as f64) - q.eval(u)).hypot();
            best = best.max(d);
        }
    }
    best * (n * n * n) as f64 / d3(c).hypot()
}

fn max_implied_dev(segs: &[Seg]) -> f64 {
    let mut worst: f64 = 0.0;
    for i in 0..segs.len() - 1 {
        worst = worst.max((segs[i].p1.midpoint(segs[i + 1].p1) - segs[i].p2).hypot());
    }
    worst
}

fn max_on_curve_dev(c: &CubicBez, segs: &[Seg]) -> f64 {
    let n = segs.len();
    let mut worst: f64 = 0.0;
    for i in 0..n - 1 {
        worst = worst.max((c.eval((i + 1) as f64 / n as f64) - segs[i].p2).hypot());
    }
    worst
}

pub fn print_consts() {
    // The error of every scheme in the family is d3 times a function of the
    // parameter only, so one generic cubic gives the constants.
    let c = CubicBez::new((0.0, 0.0), (37.0, 211.0), (301.0, 250.0), (400.0, -40.0));
    let c2 = CubicBez::new((13.0, -7.0), (-120.0, 95.0), (260.0, 310.0), (95.0, 20.0));
    println!("equal-parameter error constants: max|C - S| * n^3 / |d3|  (two unrelated cubics)");
    println!("scheme                         n   const(cubic A)  const(cubic B)  closed form");
    for n in 1..=12usize {
        let s = vec![0.0; n];
        println!(
            "mid-point (s=0)             {:4}   {:.9}     {:.9}     sqrt(3)/36 = {:.9}",
            n,
            eq_const(&c, &family_spline(&c, &s)),
            eq_const(&c2, &family_spline(&c2, &s)),
            S3_36
        );
    }
    for n in 2..=12usize {
        println!(
            "cu2qu ramp                  {:4}   {:.9}     {:.9}     Phi(n) = {:.9}  limit |d3|/tol = {:.3}",
            n,
            eq_const(&c, &ramp_spline(&c, n)),
            eq_const(&c2, &ramp_spline(&c2, n)),
            phi(n),
            (n * n * n) as f64 / phi(n)
        );
    }
    for m in 1..=6usize {
        let n = 2 * m;
        println!(
            "Truong x m pieces (+1,-1,..){:4}   {:.9}     {:.9}     4/27 = {:.9}  limit |d3|/tol = {:.3}",
            n,
            eq_const(&c, &truong_m(&c, m)),
            eq_const(&c2, &truong_m(&c2, m)),
            4.0 / 27.0,
            (n * n * n) as f64 * 27.0 / 4.0
        );
    }
    println!(
        "Triolo/Groleau fixed 4         4   {:.9}     {:.9}",
        eq_const(&c, &triolo4_fixed(&c, 0.0)),
        eq_const(&c2, &triolo4_fixed(&c2, 0.0))
    );
    for n in [2usize, 4, 8, 16] {
        let mut s = vec![0.0; n];
        s[0] = 1.0;
        s[n - 1] = -1.0;
        // Skia keeps the cubic's own points as on-curve points, so its pieces
        // are separate quadratics; the constant is the worst piece's.
        let mut worst: f64 = 0.0;
        for k in 0..n {
            let p = piece(&c, k, n);
            let ctrl = start_tangent_point(&p).lerp(end_tangent_point(&p), (1.0 - s[k]) / 2.0);
            let q = kurbo::QuadBez::new(p.p0, ctrl, p.p3);
            for j in 0..=20000 {
                let u = j as f64 / 20000.0;
                worst = worst.max((p.eval(u) - q.eval(u)).hypot());
            }
        }
        println!(
            "Skia pattern (+1,0,..,0,-1) {:4}   {:.9}                     4/27 = {:.9} (first and last piece)",
            n,
            worst * (n * n * n) as f64 / d3(&c).hypot(),
            4.0 / 27.0
        );
    }
    println!();
    println!("identities (largest deviation, font units, on cubic A):");
    let t2 = truong2_fixed(&c, 0.0);
    let r2 = ramp_spline(&c, 2);
    let dev = t2.iter().zip(r2.iter()).map(|(a, b)| (a.p1 - b.p1).hypot().max((a.p2 - b.p2).hypot())).fold(0.0, f64::max);
    println!("  Truong (gamma=1/2) vs cu2qu ramp n=2:                 {:.3e}", dev);
    let t4 = triolo4_fixed(&c, 0.0);
    let f4 = family_spline(&c, &[1.0, -1.0, 1.0, -1.0]);
    let dev = t4.iter().zip(f4.iter()).map(|(a, b)| (a.p1 - b.p1).hypot().max((a.p2 - b.p2).hypot())).fold(0.0, f64::max);
    println!("  Triolo/Groleau fixed 4 vs family s=(+1,-1,+1,-1):     {:.3e}", dev);
    let m2 = truong_m(&c, 2);
    let dev = t4.iter().zip(m2.iter()).map(|(a, b)| (a.p1 - b.p1).hypot().max((a.p2 - b.p2).hypot())).fold(0.0, f64::max);
    println!("  Triolo/Groleau fixed 4 vs Truong on two halves:       {:.3e}", dev);
    println!("  Triolo/Groleau fixed 4: on-curve implied (max dev):   {:.3e}", max_implied_dev(&t4));
    for n in [2usize, 3, 5, 8] {
        let q: Vec<Seg> = c.to_quads(1e-9).map(|(_, _, q)| Seg::quad(q.p0, q.p1, q.p2)).collect();
        let _ = q;
        let s0 = family_spline(&c, &vec![0.0; n]);
        println!(
            "  mid-point spline n={}: on-curve points on the cubic:  {:.3e}  (implied by construction)",
            n,
            max_on_curve_dev(&c, &s0)
        );
    }
    for n in [3usize, 4, 7] {
        let r = ramp_spline(&c, n);
        println!(
            "  cu2qu ramp n={}: on-curve points off the cubic by      {:.4} = {:.4} |d3|/n^3",
            n,
            max_on_curve_dev(&c, &r),
            max_on_curve_dev(&c, &r) * (n * n * n) as f64 / d3(&c).hypot()
        );
    }
}
