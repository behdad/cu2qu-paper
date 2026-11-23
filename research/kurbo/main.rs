// Bit-exact comparison harness: reads the inputs written by kdiff_py.py on
// stdin, converts each with kurbo's port of cu2qu, and writes results in the
// same format as kdiff_py.py's `run` mode, every float as the 16 hex digits
// of its bits.

use std::io::{self, BufRead, BufWriter, Write};

use kurbo::{cubics_to_quadratic_splines, CubicBez, QuadSpline};

fn float(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).unwrap())
}

fn fmt(spline: &QuadSpline) -> String {
    spline
        .points()
        .iter()
        .map(|p| format!("{:016x} {:016x}", p.x.to_bits(), p.y.to_bits()))
        .collect::<Vec<_>>()
        .join(" ")
}

fn cubic(v: &[f64]) -> CubicBez {
    CubicBez::new((v[0], v[1]), (v[2], v[3]), (v[4], v[5]), (v[6], v[7]))
}

fn main() {
    let stdin = io::stdin();
    let mut out = BufWriter::new(io::stdout().lock());
    for line in stdin.lock().lines() {
        let line = line.unwrap();
        let mut fields = line.split(' ');
        match fields.next().unwrap() {
            "S" => {
                let v: Vec<f64> = fields.map(float).collect();
                match cubic(&v[1..]).approx_spline(v[0]) {
                    Some(spline) => writeln!(out, "{}", fmt(&spline)).unwrap(),
                    None => writeln!(out, "-").unwrap(),
                }
            }
            "M" => {
                let k: usize = fields.next().unwrap().parse().unwrap();
                let v: Vec<f64> = fields.map(float).collect();
                let curves: Vec<CubicBez> = (0..k).map(|i| cubic(&v[1 + 8 * i..])).collect();
                match cubics_to_quadratic_splines(&curves, v[0]) {
                    Some(splines) => {
                        let s: Vec<String> = splines.iter().map(fmt).collect();
                        writeln!(out, "{}", s.join(" | ")).unwrap();
                    }
                    None => writeln!(out, "-").unwrap(),
                }
            }
            other => panic!("unknown record {other}"),
        }
    }
}
