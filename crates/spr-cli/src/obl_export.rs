//! Export of the final LP basis near eps* for scripts/closed_form.py.

use serde_json::json;
use spr_oblivious::lp::SymModel;
use spr_oblivious::rational::fmt;
use spr_oblivious::search::EpsResult;
use std::path::Path;

/// Writes the symmetric model restricted to what the closed-form step needs: every
/// basic column and every column in the support, with its `(row, ups, members)`
/// entries (coefficient `ups - eps * members`), the nonbasic (tight) rows of the final
/// basis, and the float solution.
pub fn export(r: &EpsResult, out: &Path) {
    let m = SymModel::build(r.space);
    // The basis that produced the certified lower bound, if the exact polish won.
    let basis = r.lower_basis.as_ref().unwrap_or(&r.basis);
    let sum_row = m.rows;
    let basic: Vec<usize> = (0..m.num_cols()).filter(|&o| basis.col_status.get(o) == Some(&1)).collect();
    let support: Vec<usize> = (0..m.num_cols()).filter(|&o| r.x[o] > 1e-12).collect();
    let mut cols: Vec<usize> = basic.iter().chain(support.iter()).cloned().collect();
    cols.sort();
    cols.dedup();
    let col_json: Vec<_> = cols
        .iter()
        .map(|&o| {
            json!({
                "orbit": o,
                "rep": m.rep(o),
                "x": r.x[o],
                "basic": basis.col_status.get(o) == Some(&1),
                "entries": m.cols[o].iter().map(|&(row, u, mem)| json!([row, u, mem])).collect::<Vec<_>>(),
            })
        })
        .collect();
    let nonbasic_rows: Vec<usize> =
        (0..=m.rows).filter(|&h| basis.row_status.get(h).is_some_and(|&st| st != 1)).collect();
    let prefixes: Vec<_> = (0..m.rows).map(|h| m.space.prefix_from_index(h)).collect();
    let v = json!({
        "n": m.space.n, "s": m.space.s,
        "rows": m.rows, "sum_row": sum_row,
        "eps_at": r.x_at, "eps_of_x": r.x_eps, "a": r.a, "b": r.b,
        "lower": fmt(&r.lower), "upper": fmt(&r.upper.value), "upper_certified": r.upper.certified,
        "basis_available": !basis.col_status.is_empty(),
        "basis_source": if r.lower_basis.is_some() { "exact polish" } else { "best float solve" },
        "num_basic_cols": basic.len(),
        "nonbasic_rows": nonbasic_rows,
        "prefixes": prefixes,
        "cols": col_json,
    });
    std::fs::write(out, serde_json::to_string(&v).unwrap()).unwrap();
    println!("exported {} columns ({} basic, {} in support), {} nonbasic rows to {}",
        cols.len(), basic.len(), support.len(), v["nonbasic_rows"].as_array().unwrap().len(), out.display());
}
