// Survey harness: cubic -> quadratic conversion, real library code plus
// labelled reimplementations, and one common evaluator.
//
//   c2q-survey corpus2bin <corpus.json> <corpus.bin>
//   c2q-survey gen <method> <tol> <corpus.json> <out.bin>
//   c2q-survey eval <label> <tol> <corpus.json> <results.bin> [percubic.out]
//   c2q-survey time <method> <tol> <corpus.json> [reps]
//   c2q-survey consts
//
// Results file: "C2QRES01", u32 count, then per cubic: u32 nseg and
// nseg * { u32 kind (1 line, 2 quad), f64 x0 y0 cx cy x1 y1 }, little endian.

use kurbo::{CubicBez, Line, ParamCurve, ParamCurveExtrema, ParamCurveNearest, Point, QuadBez, Vec2};
use std::io::{BufReader, BufWriter, Read, Write};

mod reimpl;

#[derive(Clone, Copy, Debug)]
pub struct Seg {
    pub kind: u32,
    pub p0: Point,
    pub p1: Point,
    pub p2: Point,
}

impl Seg {
    pub fn quad(p0: Point, p1: Point, p2: Point) -> Seg {
        Seg { kind: 2, p0, p1, p2 }
    }
    pub fn line(p0: Point, p2: Point) -> Seg {
        Seg { kind: 1, p0, p1: p0.midpoint(p2), p2 }
    }
    fn eval(&self, t: f64) -> Point {
        if self.kind == 1 {
            self.p0.lerp(self.p2, t)
        } else {
            QuadBez::new(self.p0, self.p1, self.p2).eval(t)
        }
    }
}

fn load_corpus(path: &str) -> Vec<CubicBez> {
    let f = std::fs::File::open(path).expect("open corpus");
    let v: Vec<[[f64; 2]; 4]> = serde_json::from_reader(BufReader::new(f)).expect("parse corpus");
    v.iter()
        .map(|c| CubicBez::new((c[0][0], c[0][1]), (c[1][0], c[1][1]), (c[2][0], c[2][1]), (c[3][0], c[3][1])))
        .collect()
}

fn write_results(path: &str, res: &[Vec<Seg>]) {
    let f = std::fs::File::create(path).expect("create results");
    let mut w = BufWriter::new(f);
    w.write_all(b"C2QRES01").unwrap();
    w.write_all(&(res.len() as u32).to_le_bytes()).unwrap();
    for segs in res {
        w.write_all(&(segs.len() as u32).to_le_bytes()).unwrap();
        for s in segs {
            w.write_all(&s.kind.to_le_bytes()).unwrap();
            for v in [s.p0.x, s.p0.y, s.p1.x, s.p1.y, s.p2.x, s.p2.y] {
                w.write_all(&v.to_le_bytes()).unwrap();
            }
        }
    }
}

fn read_results(path: &str) -> Vec<Vec<Seg>> {
    let mut buf = Vec::new();
    std::fs::File::open(path).expect("open results").read_to_end(&mut buf).unwrap();
    assert_eq!(&buf[0..8], b"C2QRES01");
    let mut pos = 8;
    let rd_u32 = |pos: &mut usize| {
        let v = u32::from_le_bytes(buf[*pos..*pos + 4].try_into().unwrap());
        *pos += 4;
        v
    };
    let rd_f64 = |pos: &mut usize| {
        let v = f64::from_le_bytes(buf[*pos..*pos + 8].try_into().unwrap());
        *pos += 8;
        v
    };
    let count = rd_u32(&mut pos) as usize;
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        let n = rd_u32(&mut pos) as usize;
        let mut segs = Vec::with_capacity(n);
        for _ in 0..n {
            let kind = rd_u32(&mut pos);
            let x0 = rd_f64(&mut pos);
            let y0 = rd_f64(&mut pos);
            let cx = rd_f64(&mut pos);
            let cy = rd_f64(&mut pos);
            let x1 = rd_f64(&mut pos);
            let y1 = rd_f64(&mut pos);
            segs.push(Seg { kind, p0: Point::new(x0, y0), p1: Point::new(cx, cy), p2: Point::new(x1, y1) });
        }
        out.push(segs);
    }
    assert_eq!(pos, buf.len());
    out
}

// ---------------------------------------------------------------- generators

fn gen_kurbo_to_quads(c: &CubicBez, tol: f64) -> Vec<Seg> {
    c.to_quads(tol).map(|(_, _, q)| Seg::quad(q.p0, q.p1, q.p2)).collect()
}

fn gen_kurbo_approx_spline(c: &CubicBez, tol: f64) -> Vec<Seg> {
    match c.approx_spline(tol) {
        Some(sp) => sp.to_quads().map(|q| Seg::quad(q.p0, q.p1, q.p2)).collect(),
        None => vec![],
    }
}

fn gen_kurbo_cubics_to_quadratic_splines(c: &CubicBez, tol: f64) -> Vec<Seg> {
    match kurbo::cubics_to_quadratic_splines(&[*c], tol) {
        Some(sps) => sps[0].to_quads().map(|q| Seg::quad(q.p0, q.p1, q.p2)).collect(),
        None => vec![],
    }
}

fn gen_lyon(c: &CubicBez, tol: f64) -> Vec<Seg> {
    use lyon_geom::{point, CubicBezierSegment};
    let seg = CubicBezierSegment {
        from: point(c.p0.x, c.p0.y),
        ctrl1: point(c.p1.x, c.p1.y),
        ctrl2: point(c.p2.x, c.p2.y),
        to: point(c.p3.x, c.p3.y),
    };
    let mut out = Vec::new();
    seg.for_each_quadratic_bezier(tol, &mut |q| {
        out.push(Seg::quad(
            Point::new(q.from.x, q.from.y),
            Point::new(q.ctrl.x, q.ctrl.y),
            Point::new(q.to.x, q.to.y),
        ));
    });
    out
}

fn gen_lyon017(c: &CubicBez, tol: f64) -> Vec<Seg> {
    use lyon_geom_017::{cubic_to_quadratic::cubic_to_quadratics, point, CubicBezierSegment};
    let seg = CubicBezierSegment {
        from: point(c.p0.x, c.p0.y),
        ctrl1: point(c.p1.x, c.p1.y),
        ctrl2: point(c.p2.x, c.p2.y),
        to: point(c.p3.x, c.p3.y),
    };
    let mut out = Vec::new();
    cubic_to_quadratics(&seg, tol, &mut |q| {
        out.push(Seg::quad(
            Point::new(q.from.x, q.from.y),
            Point::new(q.ctrl.x, q.ctrl.y),
            Point::new(q.to.x, q.to.y),
        ));
    });
    out
}

fn gen_lyon017_monotonic(c: &CubicBez, tol: f64) -> Vec<Seg> {
    use lyon_geom_017::{cubic_to_quadratic::cubic_to_monotonic_quadratics, point, CubicBezierSegment};
    let seg = CubicBezierSegment {
        from: point(c.p0.x, c.p0.y),
        ctrl1: point(c.p1.x, c.p1.y),
        ctrl2: point(c.p2.x, c.p2.y),
        to: point(c.p3.x, c.p3.y),
    };
    let mut out = Vec::new();
    cubic_to_monotonic_quadratics(&seg, tol, &mut |m| {
        let q = m.segment();
        out.push(Seg::quad(
            Point::new(q.from.x, q.from.y),
            Point::new(q.ctrl.x, q.ctrl.y),
            Point::new(q.to.x, q.to.y),
        ));
    });
    out
}

// Oracles (not implementations): the least n for which a given construction is
// within tolerance in the evaluator's own Hausdorff measure.  They show what an
// exact geometric test would buy each construction.

fn haus(c: &CubicBez, segs: &[Seg]) -> f64 {
    let r = eval_one(c, segs);
    r.haus_c2s.max(r.haus_s2c)
}

fn oracle_ramp(c: &CubicBez, tol: f64) -> Vec<Seg> {
    if let Some(s) = reimpl::tangent_single(c) {
        if haus(c, &s) <= tol {
            return s;
        }
    }
    for n in 2..=100 {
        let s = reimpl::ramp_spline(c, n);
        if haus(c, &s) <= tol {
            return s;
        }
    }
    vec![]
}

fn oracle_qm(c: &CubicBez, tol: f64) -> Vec<Seg> {
    for n in 1..=100 {
        let s = reimpl::family_spline(c, &vec![0.0; n]);
        if haus(c, &s) <= tol {
            return s;
        }
    }
    vec![]
}

fn oracle_tangent(c: &CubicBez, tol: f64) -> Vec<Seg> {
    for n in 1..=24 {
        if let Some(s) = reimpl::tangent_uniform(c, n) {
            if haus(c, &s) <= tol {
                return s;
            }
        }
    }
    oracle_ramp(c, tol)
}

fn oracle_tangent_implied(c: &CubicBez, tol: f64) -> Vec<Seg> {
    for n in 1..=24 {
        if let Some(s) = reimpl::tangent_implied(c, n) {
            if haus(c, &s) <= tol {
                return s;
            }
        }
    }
    oracle_ramp(c, tol)
}

fn oracle_truong(c: &CubicBez, tol: f64) -> Vec<Seg> {
    if let Some(s) = reimpl::tangent_single(c) {
        if haus(c, &s) <= tol {
            return s;
        }
    }
    for m in 1..=50 {
        let s = reimpl::truong_m(c, m);
        if haus(c, &s) <= tol {
            return s;
        }
    }
    vec![]
}

type GenFn = fn(&CubicBez, f64) -> Vec<Seg>;

fn method(name: &str) -> GenFn {
    match name {
        // real library code
        "kurbo_to_quads" => gen_kurbo_to_quads,
        "kurbo_approx_spline" => gen_kurbo_approx_spline,
        "kurbo_cubics_to_quadratic_splines" => gen_kurbo_cubics_to_quadratic_splines,
        "lyon" => gen_lyon,
        "lyon017" => gen_lyon017,
        "lyon017_monotonic" => gen_lyon017_monotonic,
        // reimplementations (see reimpl.rs)
        "reimpl_hb_direct" => reimpl::hb_direct,
        "reimpl_cu2qu_trial" => reimpl::cu2qu_trial,
        "reimpl_colomitchi_half" => reimpl::colomitchi_half,
        "reimpl_colomitchi_adaptive" => reimpl::colomitchi_adaptive,
        "reimpl_colomitchi_symmetric" => reimpl::colomitchi_symmetric,
        "reimpl_truong_presub" => reimpl::truong_presub,
        "reimpl_truong2_fixed" => reimpl::truong2_fixed,
        "reimpl_triolo4_fixed" => reimpl::triolo4_fixed,
        "reimpl_groleau_bspline2_fixed" => reimpl::groleau_bspline2_fixed,
        "reimpl_alfke_penner" => reimpl::alfke_penner,
        "oracle_ramp" => oracle_ramp,
        "oracle_qm" => oracle_qm,
        "oracle_truong" => oracle_truong,
        "oracle_tangent" => oracle_tangent,
        "oracle_tangent_implied" => oracle_tangent_implied,
        _ => panic!("unknown method {name}"),
    }
}

// ----------------------------------------------------------------- evaluator

fn golden_max<F: FnMut(f64) -> f64>(mut f: F, mut a: f64, mut b: f64, iters: usize) -> f64 {
    let gr = 0.618_033_988_749_894_9;
    let mut c = b - gr * (b - a);
    let mut d = a + gr * (b - a);
    let mut fc = f(c);
    let mut fd = f(d);
    let mut best = fc.max(fd);
    for _ in 0..iters {
        if fc > fd {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = f(d);
        }
        best = best.max(fc).max(fd);
    }
    best
}

fn golden_min<F: FnMut(f64) -> f64>(mut f: F, a: f64, b: f64, iters: usize) -> f64 {
    -golden_max(|t| -f(t), a, b, iters)
}

/// Distance from p to one segment: the analytic nearest point (kurbo), guarded by
/// a sampled-and-refined estimate in case the root finder misses a root.
fn dist_to_seg(s: &Seg, p: Point) -> f64 {
    if s.kind == 1 {
        return Line::new(s.p0, s.p2).nearest(p, 1e-12).distance_sq.sqrt();
    }
    let q = QuadBez::new(s.p0, s.p1, s.p2);
    let a = q.nearest(p, 1e-12).distance_sq;
    const K: usize = 8;
    let mut best = f64::INFINITY;
    let mut bi = 0;
    for i in 0..=K {
        let d = (q.eval(i as f64 / K as f64) - p).hypot2();
        if d < best {
            best = d;
            bi = i;
        }
    }
    let lo = (bi as f64 - 1.0).max(0.0) / K as f64;
    let hi = (bi as f64 + 1.0).min(K as f64) / K as f64;
    let b = golden_min(|t| (q.eval(t) - p).hypot2(), lo, hi, 40).min(best);
    a.min(b).sqrt()
}

fn dist_to_spline(segs: &[Seg], p: Point) -> f64 {
    let mut best = f64::INFINITY;
    for s in segs {
        let d = dist_to_seg(s, p);
        if d < best {
            best = d;
        }
    }
    best
}

struct CubicSampler {
    c: CubicBez,
    pts: Vec<Point>,
}

impl CubicSampler {
    fn new(c: &CubicBez, k: usize) -> Self {
        let pts = (0..=k).map(|i| c.eval(i as f64 / k as f64)).collect();
        CubicSampler { c: *c, pts }
    }
    fn dist(&self, p: Point) -> f64 {
        let k = self.pts.len() - 1;
        let mut best = f64::INFINITY;
        let mut bi = 0;
        for (i, q) in self.pts.iter().enumerate() {
            let d = (*q - p).hypot2();
            if d < best {
                best = d;
                bi = i;
            }
        }
        let lo = (bi as f64 - 1.0).max(0.0) / k as f64;
        let hi = (bi as f64 + 1.0).min(k as f64) / k as f64;
        let c = self.c;
        golden_min(|t| (c.eval(t) - p).hypot2(), lo, hi, 32).min(best).sqrt()
    }
}

fn angle_deg(a: Vec2, b: Vec2) -> f64 {
    if a.hypot2() == 0.0 || b.hypot2() == 0.0 {
        return 0.0;
    }
    a.cross(b).abs().atan2(a.dot(b)).to_degrees()
}

fn cubic_start_dir(c: &CubicBez) -> Vec2 {
    let mut d = c.p1 - c.p0;
    if d.hypot2() == 0.0 {
        d = c.p2 - c.p0;
    }
    if d.hypot2() == 0.0 {
        d = c.p3 - c.p0;
    }
    d
}

fn cubic_end_dir(c: &CubicBez) -> Vec2 {
    let mut d = c.p3 - c.p2;
    if d.hypot2() == 0.0 {
        d = c.p3 - c.p1;
    }
    if d.hypot2() == 0.0 {
        d = c.p3 - c.p0;
    }
    d
}

fn seg_start_dir(s: &Seg) -> Vec2 {
    let mut d = s.p1 - s.p0;
    if s.kind == 1 || d.hypot2() == 0.0 {
        d = s.p2 - s.p0;
    }
    d
}

fn seg_end_dir(s: &Seg) -> Vec2 {
    let mut d = s.p2 - s.p1;
    if s.kind == 1 || d.hypot2() == 0.0 {
        d = s.p2 - s.p0;
    }
    d
}

#[derive(Clone, Copy, Default, Debug)]
struct PerCubic {
    nseg: u32,
    nquad: u32,
    tt_points: u32,
    explicit_on: u32,
    haus_c2s: f64,
    haus_s2c: f64,
    eq_err: f64,
    tan_start: f64,
    tan_end: f64,
    kink: f64,
    gap: f64,
    bbox: f64,
    d3: f64,
    degenerate: bool,
}

const SAMPLES_PER_SEG: usize = 32;

fn eval_one(c: &CubicBez, segs: &[Seg]) -> PerCubic {
    let mut r = PerCubic::default();
    let n = segs.len();
    r.nseg = n as u32;
    r.nquad = segs.iter().filter(|s| s.kind == 2).count() as u32;
    r.d3 = (c.p3.to_vec2() - 3.0 * c.p2.to_vec2() + 3.0 * c.p1.to_vec2() - c.p0.to_vec2()).hypot();
    if n == 0 {
        r.haus_c2s = f64::NAN;
        r.haus_s2c = f64::NAN;
        r.eq_err = f64::NAN;
        return r;
    }
    r.degenerate = c.p0 == c.p1 && c.p1 == c.p2 && c.p2 == c.p3;

    // TrueType points: off-curve points plus the interior on-curve points that are
    // not implied (not the midpoint of the off-curve points on both sides).
    let mut explicit_on = 0;
    let mut gap: f64 = (segs[0].p0 - c.p0).hypot().max((segs[n - 1].p2 - c.p3).hypot());
    let mut kink: f64 = 0.0;
    for i in 0..n - 1 {
        let a = &segs[i];
        let b = &segs[i + 1];
        gap = gap.max((a.p2 - b.p0).hypot());
        let implied = a.kind == 2 && b.kind == 2 && (a.p1.midpoint(b.p1) - a.p2).hypot() <= 1e-6;
        if !implied {
            explicit_on += 1;
        }
        kink = kink.max(angle_deg(seg_end_dir(a), seg_start_dir(b)));
    }
    r.explicit_on = explicit_on;
    r.tt_points = r.nquad + explicit_on;
    r.gap = gap;
    r.kink = kink;

    // End tangents.
    r.tan_start = angle_deg(cubic_start_dir(c), seg_start_dir(&segs[0]));
    r.tan_end = angle_deg(cubic_end_dir(c), seg_end_dir(&segs[n - 1]));

    // Bounding box.
    let bc = ParamCurveExtrema::bounding_box(c);
    let mut x0 = f64::INFINITY;
    let mut y0 = f64::INFINITY;
    let mut x1 = f64::NEG_INFINITY;
    let mut y1 = f64::NEG_INFINITY;
    for s in segs {
        let b = if s.kind == 1 {
            kurbo::Rect::from_points(s.p0, s.p2)
        } else {
            ParamCurveExtrema::bounding_box(&QuadBez::new(s.p0, s.p1, s.p2))
        };
        x0 = x0.min(b.x0);
        y0 = y0.min(b.y0);
        x1 = x1.max(b.x1);
        y1 = y1.max(b.y1);
    }
    r.bbox = (bc.x0 - x0).abs().max((bc.y0 - y0).abs()).max((bc.x1 - x1).abs()).max((bc.y1 - y1).abs());

    // Hausdorff distance, cubic -> spline.
    let m = SAMPLES_PER_SEG * n.max(2);
    let vals: Vec<f64> = (0..=m).map(|j| dist_to_spline(segs, c.eval(j as f64 / m as f64))).collect();
    r.haus_c2s = refine_peaks(&vals, |t| dist_to_spline(segs, c.eval(t)), 20);

    // Hausdorff distance, spline -> cubic.
    let sampler = CubicSampler::new(c, m.max(64));
    let mut worst: f64 = 0.0;
    let per_seg: Vec<Vec<f64>> = segs
        .iter()
        .map(|s| (0..=SAMPLES_PER_SEG).map(|j| sampler.dist(s.eval(j as f64 / SAMPLES_PER_SEG as f64))).collect())
        .collect();
    let overall = per_seg.iter().flat_map(|v| v.iter().copied()).fold(0.0, f64::max);
    for (s, vals) in segs.iter().zip(per_seg.iter()) {
        worst = worst.max(refine_peaks(vals, |u| sampler.dist(s.eval(u)), 20));
    }
    r.haus_s2c = worst.max(overall);

    // Equal-parameter distance, taking segment k to cover [k/n, (k+1)/n] of the cubic.
    // Only meaningful for the methods that split uniformly.
    let f = |k: usize, u: f64| (c.eval((k as f64 + u) / n as f64) - segs[k].eval(u)).hypot();
    let mut worst: f64 = 0.0;
    for k in 0..n {
        let vals: Vec<f64> = (0..=64).map(|j| f(k, j as f64 / 64.0)).collect();
        worst = worst.max(refine_peaks(&vals, |u| f(k, u), 30));
    }
    r.eq_err = worst;
    r
}

/// The maximum of a function sampled at vals[j] = f(j / m): every local peak of
/// the samples within 10% of the largest is refined by golden-section search,
/// and so are the first and the last interval whatever their samples say (a
/// segment whose control point lies behind its end point overshoots that end
/// in a hook much narrower than the sampling step).
fn refine_peaks<F: FnMut(f64) -> f64>(vals: &[f64], mut f: F, iters: usize) -> f64 {
    let m = vals.len() - 1;
    let best = vals.iter().copied().fold(0.0, f64::max);
    let mut out = best;
    for j in 0..=m {
        let v = vals[j];
        if v < 0.9 * best {
            continue;
        }
        let left = if j > 0 { vals[j - 1] } else { f64::NEG_INFINITY };
        let right = if j < m { vals[j + 1] } else { f64::NEG_INFINITY };
        if v >= left && v >= right {
            let lo = (j as f64 - 1.0).max(0.0) / m as f64;
            let hi = (j as f64 + 1.0).min(m as f64) / m as f64;
            out = out.max(golden_max(&mut f, lo, hi, iters));
        }
    }
    // the two end intervals, at geometrically shrinking scales
    let h = 1.0 / m as f64;
    let mut w = h;
    for _ in 0..6 {
        out = out.max(golden_max(&mut f, 0.0, w, 12));
        out = out.max(golden_max(&mut f, 1.0 - w, 1.0, 12));
        w *= 0.125;
    }
    out
}

/// Independent check of the Hausdorff distance: both curves are sampled densely
/// and every sample's distance to the other curve's polyline is taken by brute
/// force.  No root finding, no refinement.
fn brute_hausdorff(c: &CubicBez, segs: &[Seg], per_seg: usize) -> (f64, f64) {
    let n = segs.len().max(1);
    let mc = per_seg * n;
    let cp: Vec<Point> = (0..=mc).map(|j| c.eval(j as f64 / mc as f64)).collect();
    let mut sp: Vec<Point> = Vec::with_capacity(per_seg * n + 1);
    for (k, s) in segs.iter().enumerate() {
        let start = if k == 0 { 0 } else { 1 };
        for j in start..=per_seg {
            sp.push(s.eval(j as f64 / per_seg as f64));
        }
    }
    fn one_sided(a: &[Point], b: &[Point]) -> f64 {
        let mut worst: f64 = 0.0;
        for p in a {
            let mut best = f64::INFINITY;
            for w in b.windows(2) {
                let ab = w[1] - w[0];
                let l2 = ab.hypot2();
                let u = if l2 > 0.0 { ((*p - w[0]).dot(ab) / l2).clamp(0.0, 1.0) } else { 0.0 };
                let d = (*p - (w[0] + ab * u)).hypot2();
                if d < best {
                    best = d;
                }
            }
            worst = worst.max(best);
        }
        worst.sqrt()
    }
    (one_sided(&cp, &sp), one_sided(&sp, &cp))
}

/// Same n for every construction: the Hausdorff and equal-parameter distances of
/// each, as multiples of |d3| / n^3, over the cubics that are not close to
/// quadratic already (|d3| >= 8).
fn cmd_samen(corpus: &str) {
    let cubics = load_corpus(corpus);
    let d3 = |c: &CubicBez| (c.p3.to_vec2() - 3.0 * c.p2.to_vec2() + 3.0 * c.p1.to_vec2() - c.p0.to_vec2()).hypot();
    let idx: Vec<usize> = (0..cubics.len()).filter(|&i| d3(&cubics[i]) >= 8.0).collect();
    println!("{} of {} cubics have |d3| >= 8", idx.len(), cubics.len());
    println!(
        "{:>2} {:28} {:>8} {:>10} {:>10} {:>10} {:>10} {:>10} {:>12}",
        "n", "construction", "cubics", "eq const", "H med", "H mean", "H p90", "H max", "H/H(mid-pt)"
    );
    for n in 1..=6usize {
        let names: Vec<&str> = if n == 1 {
            vec!["mid-point (s=0)", "tangent intersection"]
        } else if n % 2 == 0 {
            vec![
                "mid-point (s=0)",
                "cu2qu ramp",
                "Truong pairs (+1,-1,..)",
                "tangent intersection",
                "tangent int., implied on",
                "Skia pattern (+1,0..0,-1)",
            ]
        } else {
            vec!["mid-point (s=0)", "cu2qu ramp", "tangent intersection", "tangent int., implied on", "Skia pattern (+1,0..0,-1)"]
        };
        let rows = par_map(idx.len(), |k| {
            let c = &cubics[idx[k]];
            let scale = (n * n * n) as f64 / d3(c);
            let mut out: Vec<Option<(f64, f64)>> = Vec::new();
            for name in &names {
                let segs: Option<Vec<Seg>> = match *name {
                    "mid-point (s=0)" => Some(reimpl::family_spline(c, &vec![0.0; n])),
                    "cu2qu ramp" => Some(reimpl::ramp_spline(c, n)),
                    "Truong pairs (+1,-1,..)" => Some(reimpl::truong_m(c, n / 2)),
                    "tangent intersection" => reimpl::tangent_uniform(c, n),
                    "tangent int., implied on" => reimpl::tangent_implied(c, n),
                    "Skia pattern (+1,0..0,-1)" => Some(reimpl::skia_pattern(c, n)),
                    _ => unreachable!(),
                };
                out.push(segs.map(|s| {
                    let e = eval_one(c, &s);
                    (e.haus_c2s.max(e.haus_s2c) * scale, e.eq_err * scale)
                }));
            }
            out
        });
        for (j, name) in names.iter().enumerate() {
            let mut h: Vec<f64> = Vec::new();
            let mut eq: Vec<f64> = Vec::new();
            let mut ratio: Vec<f64> = Vec::new();
            for r in &rows {
                if let Some((hh, ee)) = r[j] {
                    h.push(hh);
                    eq.push(ee);
                    if let Some((h0, _)) = r[0] {
                        if h0 > 1e-9 {
                            ratio.push(hh / h0);
                        }
                    }
                }
            }
            h.sort_by(|a, b| a.partial_cmp(b).unwrap());
            eq.sort_by(|a, b| a.partial_cmp(b).unwrap());
            ratio.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mean = h.iter().sum::<f64>() / h.len() as f64;
            println!(
                "{:>2} {:28} {:>8} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>12.3}",
                n,
                name,
                h.len(),
                pct(&eq, 0.5),
                pct(&h, 0.5),
                mean,
                pct(&h, 0.9),
                h.last().copied().unwrap_or(f64::NAN),
                pct(&ratio, 0.5)
            );
        }
    }
}

/// Convergence order in the Hausdorff measure: error at n and 2n for a few
/// inflection-free cubics.
fn cmd_rate() {
    let cs = [
        CubicBez::new((0.0, 0.0), (0.0, 55.0), (45.0, 100.0), (100.0, 100.0)),
        CubicBez::new((300.0, 0.0), (447.0, 0.0), (547.0, 63.0), (547.0, 192.0)),
        CubicBez::new((0.0, 0.0), (37.0, 211.0), (301.0, 250.0), (400.0, -40.0)),
    ];
    for c in &cs {
        println!("cubic {:?} {:?} {:?} {:?}", c.p0, c.p1, c.p2, c.p3);
        println!("{:>4} {:>14} {:>7} {:>14} {:>7} {:>14} {:>7}", "n", "mid-point H", "order", "tangent-int H", "order", "cu2qu ramp H", "order");
        let mut prev: Option<(f64, f64, f64)> = None;
        for n in [2usize, 4, 8, 16, 32, 64] {
            let a = haus(c, &reimpl::family_spline(c, &vec![0.0; n]));
            let b = reimpl::tangent_uniform(c, n).map(|s| haus(c, &s)).unwrap_or(f64::NAN);
            let r = haus(c, &reimpl::ramp_spline(c, n));
            let ord = |x: f64, p: f64| (p / x).log2();
            match prev {
                Some((pa, pb, pr)) => println!(
                    "{:>4} {:>14.6e} {:>7.2} {:>14.6e} {:>7.2} {:>14.6e} {:>7.2}",
                    n, a, ord(a, pa), b, ord(b, pb), r, ord(r, pr)
                ),
                None => println!("{:>4} {:>14.6e} {:>7} {:>14.6e} {:>7} {:>14.6e} {:>7}", n, a, "", b, "", r, ""),
            }
            prev = Some((a, b, r));
        }
    }
}

fn cmd_verify(corpus: &str, results: &str, count: usize) {
    let cubics = load_corpus(corpus);
    let res = read_results(results);
    // a fixed pseudo-random subset
    let mut idxs = Vec::new();
    let mut x: u64 = 0x9E3779B97F4A7C15;
    while idxs.len() < count {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        let i = (x % cubics.len() as u64) as usize;
        if !res[i].is_empty() {
            idxs.push(i);
        }
    }
    let out = par_map(idxs.len(), |k| {
        let i = idxs[k];
        let e = eval_one(&cubics[i], &res[i]);
        let (b1, b2) = brute_hausdorff(&cubics[i], &res[i], 2000);
        (e.haus_c2s, e.haus_s2c, b1, b2)
    });
    let mut worst_abs: f64 = 0.0;
    let mut worst_rel: f64 = 0.0;
    let mut worst_under: f64 = 0.0;
    let mut worst_k = 0;
    for (k, (h1, h2, b1, b2)) in out.iter().enumerate() {
        for (h, b) in [(h1, b1), (h2, b2)] {
            if (h - b).abs() > worst_abs {
                worst_abs = (h - b).abs();
                worst_k = k;
            }
            if h.max(*b) > 1e-6 {
                worst_rel = worst_rel.max((h - b).abs() / h.max(*b));
            }
            // the harness reporting less than brute force is the direction that matters
            worst_under = worst_under.max(b - h);
        }
    }
    println!(
        "{}: {} cubics: largest |harness - brute| = {:.3e} units ({:.3e} relative); harness below brute by at most {:.3e} units",
        results, count, worst_abs, worst_rel, worst_under
    );
    if std::env::var("VERIFY_VERBOSE").is_ok() {
        let i = idxs[worst_k];
        let c = &cubics[i];
        println!("  worst: cubic {} {:?} {:?} {:?} {:?}", i, c.p0, c.p1, c.p2, c.p3);
        println!("  harness c2s {:.6e} s2c {:.6e}; brute c2s {:.6e} s2c {:.6e}", out[worst_k].0, out[worst_k].1, out[worst_k].2, out[worst_k].3);
        for s in &res[i] {
            println!("  seg kind {} {:?} {:?} {:?}", s.kind, s.p0, s.p1, s.p2);
        }
    }
}

fn par_map<T: Send, F: Fn(usize) -> T + Sync>(n: usize, f: F) -> Vec<T> {
    let threads = std::thread::available_parallelism().map(|x| x.get()).unwrap_or(4).min(24);
    let chunk = (n + threads - 1) / threads;
    let mut out: Vec<Option<T>> = (0..n).map(|_| None).collect();
    std::thread::scope(|s| {
        for (ci, slice) in out.chunks_mut(chunk).enumerate() {
            let f = &f;
            s.spawn(move || {
                for (i, slot) in slice.iter_mut().enumerate() {
                    *slot = Some(f(ci * chunk + i));
                }
            });
        }
    });
    out.into_iter().map(|x| x.unwrap()).collect()
}

fn pct(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let idx = ((sorted.len() - 1) as f64 * p).round() as usize;
    sorted[idx]
}

fn stats(mut v: Vec<f64>) -> String {
    v.retain(|x| x.is_finite());
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mean = v.iter().sum::<f64>() / v.len().max(1) as f64;
    format!(
        "{{\"mean\":{:.6},\"p50\":{:.6},\"p90\":{:.6},\"p99\":{:.6},\"p999\":{:.6},\"max\":{:.6}}}",
        mean,
        pct(&v, 0.5),
        pct(&v, 0.9),
        pct(&v, 0.99),
        pct(&v, 0.999),
        v.last().copied().unwrap_or(f64::NAN)
    )
}

fn cmd_eval(label: &str, tol: f64, corpus: &str, results: &str, percubic: Option<&str>) {
    let cubics = load_corpus(corpus);
    let res = read_results(results);
    assert_eq!(cubics.len(), res.len());
    let per = par_map(cubics.len(), |i| eval_one(&cubics[i], &res[i]));

    let ncub = per.len();
    let failed = per.iter().filter(|p| p.nseg == 0).count();
    let ok: Vec<&PerCubic> = per.iter().filter(|p| p.nseg > 0).collect();
    let nseg_total: u64 = ok.iter().map(|p| p.nseg as u64).sum();
    let nquad_total: u64 = ok.iter().map(|p| p.nquad as u64).sum();
    let tt_total: u64 = ok.iter().map(|p| p.tt_points as u64).sum();
    let expl_total: u64 = ok.iter().map(|p| p.explicit_on as u64).sum();
    let mut hist = std::collections::BTreeMap::new();
    for p in &ok {
        *hist.entry(p.nseg).or_insert(0u32) += 1;
    }
    let hist_s: Vec<String> = hist.iter().map(|(k, v)| format!("\"{}\":{}", k, v)).collect();
    let haus: Vec<f64> = ok.iter().map(|p| p.haus_c2s.max(p.haus_s2c) / tol).collect();
    let viol = haus.iter().filter(|&&h| h > 1.0 + 1e-9).count();
    let viol1 = haus.iter().filter(|&&h| h > 1.01).count();
    let viol10 = haus.iter().filter(|&&h| h > 1.10).count();
    let c2s: Vec<f64> = ok.iter().map(|p| p.haus_c2s / tol).collect();
    let s2c: Vec<f64> = ok.iter().map(|p| p.haus_s2c / tol).collect();
    let eq: Vec<f64> = ok.iter().map(|p| p.eq_err / tol).collect();
    let eq_viol = eq.iter().filter(|&&h| h > 1.0 + 1e-9).count();
    // tangents: skip fully degenerate cubics
    let tan: Vec<f64> = ok.iter().filter(|p| !p.degenerate).map(|p| p.tan_start.max(p.tan_end)).collect();
    let tan_001 = tan.iter().filter(|&&a| a > 0.01).count();
    let tan_01 = tan.iter().filter(|&&a| a > 0.1).count();
    let tan_1 = tan.iter().filter(|&&a| a > 1.0).count();
    let multi: Vec<&&PerCubic> = ok.iter().filter(|p| p.nseg > 1).collect();
    let kink: Vec<f64> = multi.iter().map(|p| p.kink).collect();
    let kink_001 = kink.iter().filter(|&&a| a > 0.01).count();
    let kink_1 = kink.iter().filter(|&&a| a > 1.0).count();
    let gap = ok.iter().map(|p| p.gap).fold(0.0, f64::max);
    let bbox: Vec<f64> = ok.iter().map(|p| p.bbox).collect();
    // of the cubics converted to more than one segment
    let multi_n = multi.len();

    println!(
        "{{\"label\":\"{}\",\"tol\":{},\"cubics\":{},\"failed\":{},\"segs\":{},\"quads\":{},\"lines\":{},\"mean_n\":{:.4},\"hist\":{{{}}},\"tt_points\":{},\"tt_mean\":{:.4},\"explicit_on\":{},\"haus\":{},\"haus_c2s\":{},\"haus_s2c\":{},\"viol\":{},\"viol_1pct\":{},\"viol_10pct\":{},\"eq\":{},\"eq_viol\":{},\"tan_deg\":{},\"tan_gt_0.01\":{},\"tan_gt_0.1\":{},\"tan_gt_1\":{},\"multi\":{},\"kink_deg\":{},\"kink_gt_0.01\":{},\"kink_gt_1\":{},\"gap_max\":{:.3e},\"bbox\":{}}}",
        label,
        tol,
        ncub,
        failed,
        nseg_total,
        nquad_total,
        nseg_total - nquad_total,
        nseg_total as f64 / ok.len().max(1) as f64,
        hist_s.join(","),
        tt_total,
        tt_total as f64 / ok.len().max(1) as f64,
        expl_total,
        stats(haus),
        stats(c2s),
        stats(s2c),
        viol,
        viol1,
        viol10,
        stats(eq),
        eq_viol,
        stats(tan),
        tan_001,
        tan_01,
        tan_1,
        multi_n,
        stats(kink),
        kink_001,
        kink_1,
        gap,
        stats(bbox)
    );

    if let Some(path) = percubic {
        let f = std::fs::File::create(path).unwrap();
        let mut w = BufWriter::new(f);
        writeln!(w, "idx nseg nquad tt haus_c2s haus_s2c eq tan_start tan_end kink bbox d3").unwrap();
        for (i, p) in per.iter().enumerate() {
            writeln!(
                w,
                "{} {} {} {} {:.9} {:.9} {:.9} {:.6} {:.6} {:.6} {:.6} {:.6}",
                i, p.nseg, p.nquad, p.tt_points, p.haus_c2s, p.haus_s2c, p.eq_err, p.tan_start, p.tan_end, p.kink, p.bbox, p.d3
            )
            .unwrap();
        }
    }
}

fn cmd_gen(name: &str, tol: f64, corpus: &str, out: &str) {
    let cubics = load_corpus(corpus);
    let f = method(name);
    let res: Vec<Vec<Seg>> = if name.starts_with("oracle_") {
        par_map(cubics.len(), |i| f(&cubics[i], tol))
    } else {
        cubics.iter().map(|c| f(c, tol)).collect()
    };
    write_results(out, &res);
}

fn cmd_time(name: &str, tol: f64, corpus: &str, reps: usize) {
    let cubics = load_corpus(corpus);
    let mut best = f64::INFINITY;
    let mut check = 0.0;
    let mut nq = 0usize;
    // Library calls are timed through their own iteration interface where it
    // avoids allocating a Vec, since that is how a renderer would call them.
    for _ in 0..reps {
        let t0 = std::time::Instant::now();
        let mut acc = 0.0;
        let mut count = 0usize;
        match name {
            "kurbo_to_quads" => {
                for c in &cubics {
                    for (_, _, q) in c.to_quads(tol) {
                        acc += q.p1.x + q.p2.y;
                        count += 1;
                    }
                }
            }
            "kurbo_approx_spline" => {
                for c in &cubics {
                    if let Some(sp) = c.approx_spline(tol) {
                        let pts = sp.points();
                        acc += pts[1].x + pts[pts.len() - 1].y;
                        count += pts.len() - 2;
                    }
                }
            }
            "lyon" => {
                use lyon_geom::{point, CubicBezierSegment};
                for c in &cubics {
                    let seg = CubicBezierSegment {
                        from: point(c.p0.x, c.p0.y),
                        ctrl1: point(c.p1.x, c.p1.y),
                        ctrl2: point(c.p2.x, c.p2.y),
                        to: point(c.p3.x, c.p3.y),
                    };
                    seg.for_each_quadratic_bezier(tol, &mut |q| {
                        acc += q.ctrl.x + q.to.y;
                        count += 1;
                    });
                }
            }
            "lyon017" => {
                use lyon_geom_017::{cubic_to_quadratic::cubic_to_quadratics, point, CubicBezierSegment};
                for c in &cubics {
                    let seg = CubicBezierSegment {
                        from: point(c.p0.x, c.p0.y),
                        ctrl1: point(c.p1.x, c.p1.y),
                        ctrl2: point(c.p2.x, c.p2.y),
                        to: point(c.p3.x, c.p3.y),
                    };
                    cubic_to_quadratics(&seg, tol, &mut |q| {
                        acc += q.ctrl.x + q.to.y;
                        count += 1;
                    });
                }
            }
            _ => {
                let f = method(name);
                for c in &cubics {
                    let v = f(c, tol);
                    if let Some(s) = v.last() {
                        acc += s.p1.x + s.p2.y;
                    }
                    count += v.len();
                }
            }
        }
        let dt = t0.elapsed().as_secs_f64();
        best = best.min(dt);
        check = acc;
        nq = count;
    }
    println!(
        "{{\"method\":\"{}\",\"tol\":{},\"ns_per_cubic\":{:.1},\"quads\":{},\"check\":{:.3}}}",
        name,
        tol,
        best * 1e9 / cubics.len() as f64,
        nq,
        check
    );
}

fn cmd_corpus2bin(corpus: &str, out: &str) {
    let cubics = load_corpus(corpus);
    let f = std::fs::File::create(out).unwrap();
    let mut w = BufWriter::new(f);
    w.write_all(b"C2QCOR01").unwrap();
    w.write_all(&(cubics.len() as u32).to_le_bytes()).unwrap();
    for c in &cubics {
        for v in [c.p0.x, c.p0.y, c.p1.x, c.p1.y, c.p2.x, c.p2.y, c.p3.x, c.p3.y] {
            w.write_all(&v.to_le_bytes()).unwrap();
        }
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    match a[1].as_str() {
        "corpus2bin" => cmd_corpus2bin(&a[2], &a[3]),
        "gen" => cmd_gen(&a[2], a[3].parse().unwrap(), &a[4], &a[5]),
        "eval" => cmd_eval(&a[2], a[3].parse().unwrap(), &a[4], &a[5], a.get(6).map(|s| s.as_str())),
        "time" => cmd_time(&a[2], a[3].parse().unwrap(), &a[4], a.get(5).map(|s| s.parse().unwrap()).unwrap_or(5)),
        "consts" => reimpl::print_consts(),
        "samen" => cmd_samen(&a[2]),
        "rate" => cmd_rate(),
        "verify" => cmd_verify(&a[2], &a[3], a.get(4).map(|s| s.parse().unwrap()).unwrap_or(2000)),
        _ => panic!("unknown command"),
    }
}
