// Workload for `perf stat -e instructions:u`: read the single-cubic inputs
// (kdiff_py.py format) from a file, then convert them all, `reps` times.
// Usage: kbench <inputs> <reps>

use std::fs;

use kurbo::CubicBez;

fn float(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).unwrap())
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().unwrap();
    let reps: usize = args.next().unwrap().parse().unwrap();
    let text = fs::read_to_string(path).unwrap();
    let mut cases = Vec::new();
    for line in text.lines() {
        let mut fields = line.split(' ');
        if fields.next() != Some("S") {
            continue;
        }
        let v: Vec<f64> = fields.map(float).collect();
        cases.push((
            CubicBez::new((v[1], v[2]), (v[3], v[4]), (v[5], v[6]), (v[7], v[8])),
            v[0],
        ));
    }
    let mut points = 0usize;
    for _ in 0..reps {
        for (cubic, accuracy) in &cases {
            if let Some(spline) = std::hint::black_box(cubic).approx_spline(*accuracy) {
                points += spline.points().len();
            }
        }
    }
    println!("{} cubics, {} reps, {} points", cases.len(), reps, points);
}
