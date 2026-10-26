//! Diagnostic: where does the certified lower bound lose precision?
use spr_oblivious::exact::verify;
use spr_oblivious::search::solve;
use spr_oblivious::tuples::Space;
fn main() {
    let a: Vec<usize> = std::env::args().skip(1).map(|x| x.parse().unwrap()).collect();
    let r = solve(Space::new(a[0], a[1]), 1e-11).unwrap();
    let rep = verify(&r.lower_strategy).unwrap();
    println!("LP bracket [{:.12}, {:.12}], float eps(x) = {:.12} at eps {:.12}", r.a, r.b, r.x_eps, r.x_at);
    println!("exact lower {:.12}, worst prefix {:?}", spr_oblivious::rational::to_f64(&r.lower), rep.worst);
    let mut xs: Vec<f64> = r.x.iter().cloned().filter(|&v| v > 0.0).collect();
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    println!("smallest positive orbit masses: {:?}", &xs[..xs.len().min(8)]);
    let basic_cols = r.basis.col_status.iter().filter(|&&c| c == 1).count();
    let nonbasic_rows = r.basis.row_status.iter().filter(|&&c| c != 1).count();
    println!("basis: {basic_cols} basic structural columns, {nonbasic_rows} nonbasic rows (incl. sum row), {} rows total", r.basis.row_status.len());
}
