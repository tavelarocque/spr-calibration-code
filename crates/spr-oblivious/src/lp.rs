//! Float LPs (HiGHS) over reflection-symmetric strategies.
//!
//! Columns are reflection orbits `{t, reflect(t)}` of tuples (for `s >= 2` every orbit
//! has two members). The column variable `x_o` is the probability of EACH member, so
//! the normalisation row is `sum_o 2 x_o = 1`. For a symmetric pi, `D(h) = U(reflect h)`,
//! so the feasibility LP needs only the rows `U(h) - eps mu(h) >= 0`, one per prefix.
//!
//! Every column is stored as a list of `(row, ups, members)`: `members` of the orbit's
//! tuples extend the row's prefix and `ups` of them have every later cell above, so the
//! coefficient is `ups - eps * members`.
//!
//! Nothing produced here is trusted: every number that reaches data/ is re-checked in
//! [`crate::exact`].

use crate::tuples::{Space, up_down};
use highs::{ColProblem, HighsModelStatus, RowProblem, Sense};
use std::collections::HashMap;

pub struct SymModel {
    pub space: Space,
    /// Orbit representatives (the lexicographically smaller member), stride `s`.
    pub reps: Vec<u8>,
    /// Per column: `(row, ups, members)`.
    pub cols: Vec<Vec<(u32, u8, u8)>>,
    pub rows: usize,
}

impl SymModel {
    pub fn build(space: Space) -> Self {
        assert!(space.s >= 2);
        let all = space.all_tuples();
        let mut reps = Vec::new();
        let mut cols = Vec::new();
        for t in all.chunks(space.s) {
            let rt = space.reflect(t);
            if t > &rt[..] {
                continue;
            }
            let mut entries: Vec<(u32, u8, u8)> = Vec::with_capacity(2 * (space.s - 1));
            for member in [t, &rt[..]] {
                let (up, _) = up_down(member);
                for i in 1..space.s {
                    let row = space.prefix_index(&member[..i]) as u32;
                    let u = up[i - 1] as u8;
                    match entries.iter_mut().find(|e| e.0 == row) {
                        Some(e) => {
                            e.1 += u;
                            e.2 += 1;
                        }
                        None => entries.push((row, u, 1)),
                    }
                }
            }
            entries.sort_unstable();
            reps.extend_from_slice(t);
            cols.push(entries);
        }
        SymModel { space, reps, cols, rows: space.num_prefixes() }
    }

    pub fn num_cols(&self) -> usize {
        self.cols.len()
    }

    pub fn nnz(&self) -> usize {
        self.cols.iter().map(|c| c.len()).sum::<usize>() + self.cols.len()
    }

    pub fn rep(&self, o: usize) -> &[u8] {
        &self.reps[o * self.space.s..(o + 1) * self.space.s]
    }

    /// Float eps(pi) of the symmetric strategy `x` (probability of each orbit member),
    /// evaluated on the full U and D rows. Exact checks live in `crate::exact`.
    pub fn eps_of(&self, x: &[f64]) -> f64 {
        let mut mu = vec![0f64; self.rows];
        let mut u = vec![0f64; self.rows];
        for (o, col) in self.cols.iter().enumerate() {
            for &(row, ups, members) in col {
                mu[row as usize] += members as f64 * x[o];
                u[row as usize] += ups as f64 * x[o];
            }
        }
        let tol = 1e-14;
        let mut eps = 1f64;
        for h in 0..self.rows {
            if mu[h] > tol {
                let p = self.space.prefix_from_index(h);
                let rh = self.space.prefix_index(&self.space.reflect(&p));
                // D(h) = U(reflect h) for a symmetric strategy.
                eps = eps.min(u[h] / mu[h]).min(u[rh] / mu[h]);
            }
        }
        eps
    }
}

pub enum Feas {
    Feasible { x: Vec<f64>, basis: Basis },
    Infeasible,
    Unknown(String),
}

/// HiGHS basis statuses: 0 lower, 1 basic, 2 upper, 3 zero, 4 nonbasic.
#[derive(Clone, Debug, Default)]
pub struct Basis {
    pub col_status: Vec<i32>,
    pub row_status: Vec<i32>,
}

fn configure(model: &mut highs::Model) {
    model.make_quiet();
    model.set_option("primal_feasibility_tolerance", 1e-10);
    model.set_option("dual_feasibility_tolerance", 1e-10);
    model.set_option("threads", 1);
}

// HighsInt is i32 in the default build but may be i64; the cast keeps both working.
#[allow(clippy::unnecessary_cast)]
fn basis_of(solved: &highs::SolvedModel, ncols: usize, nrows: usize) -> Basis {
    let mut col_status = vec![0 as highs_sys::HighsInt; ncols];
    let mut row_status = vec![0 as highs_sys::HighsInt; nrows];
    // SAFETY: the buffers have exactly the model's column and row counts.
    let st = unsafe { highs_sys::Highs_getBasis(solved.as_ptr(), col_status.as_mut_ptr(), row_status.as_mut_ptr()) };
    if st != 0 && st != 1 {
        return Basis::default();
    }
    Basis {
        col_status: col_status.into_iter().map(|v| v as i32).collect(),
        row_status: row_status.into_iter().map(|v| v as i32).collect(),
    }
}

/// Is there a symmetric pi with `U(h) >= eps mu(h)` for every prefix? Near the
/// threshold HiGHS sometimes ends with status Unknown; the solve is then retried with
/// presolve off, then with the interior-point solver.
pub fn feasible(m: &SymModel, eps: f64) -> Feas {
    let mut last = Feas::Unknown(String::new());
    for variant in 0..3 {
        last = feasible_with(m, eps, variant);
        if !matches!(last, Feas::Unknown(_)) {
            return last;
        }
    }
    last
}

fn feasible_with(m: &SymModel, eps: f64, variant: u32) -> Feas {
    let mut pb = ColProblem::default();
    let rows: Vec<highs::Row> = (0..m.rows).map(|_| pb.add_row(0.0..)).collect();
    let sum = pb.add_row(1.0..=1.0);
    for col in &m.cols {
        let factors = col
            .iter()
            .map(|&(r, ups, mem)| (rows[r as usize], ups as f64 - eps * mem as f64))
            .chain(std::iter::once((sum, 2.0)));
        pb.add_column(0.0, 0.0.., factors);
    }
    let mut model = pb.optimise(Sense::Maximise);
    configure(&mut model);
    match variant {
        1 => model.set_option("presolve", "off"),
        2 => {
            model.set_option("presolve", "off");
            model.set_option("solver", "ipm");
        }
        _ => {}
    }
    let solved = match model.try_solve() {
        Ok(s) => s,
        Err(e) => return Feas::Unknown(format!("{e:?}")),
    };
    match solved.status() {
        HighsModelStatus::Optimal => {
            let x = solved.get_solution().columns().to_vec();
            let basis = basis_of(&solved, m.num_cols(), m.rows + 1);
            Feas::Feasible { x, basis }
        }
        HighsModelStatus::Infeasible => Feas::Infeasible,
        other => Feas::Unknown(format!("{other:?}")),
    }
}

/// Farkas LP for the symmetric system at `eps`: find `y >= 0` (one entry per U-row) with
/// `sum_h y_h (ups - eps members) <= -1` in every column. Minimises `sum y`.
pub fn farkas(m: &SymModel, eps: f64) -> Option<Vec<f64>> {
    let mut pb = RowProblem::default();
    let ys: Vec<highs::Col> = (0..m.rows).map(|_| pb.add_column(1.0, 0.0..)).collect();
    for col in &m.cols {
        let factors: Vec<(highs::Col, f64)> =
            col.iter().map(|&(r, ups, mem)| (ys[r as usize], ups as f64 - eps * mem as f64)).collect();
        pb.add_row(..=-1.0, factors);
    }
    let mut model = pb.optimise(Sense::Minimise);
    configure(&mut model);
    let solved = model.try_solve().ok()?;
    (solved.status() == HighsModelStatus::Optimal).then(|| solved.get_solution().columns().to_vec())
}

pub struct OvSolution {
    pub value: f64,
    pub x: Vec<f64>,
    /// Labeler probability of '+' at each prefix, from the duals.
    pub a: Vec<f64>,
}

/// ov* over symmetric strategies: maximise `1 + sum_h z_h` subject to `z_h <= U(h)` and
/// `z_h <= D(h) = U(reflect h)` for every prefix, `sum 2 x_o = 1`.
pub fn ov(m: &SymModel) -> Option<OvSolution> {
    let sp = m.space;
    let refl: Vec<usize> = (0..m.rows).map(|h| sp.prefix_index(&sp.reflect(&sp.prefix_from_index(h)))).collect();
    let mut pb = ColProblem::default();
    // Row 2h: z_h - U(h) <= 0. Row 2h + 1: z_h - U(reflect h) <= 0.
    let rows: Vec<highs::Row> = (0..2 * m.rows).map(|_| pb.add_row(..=0.0)).collect();
    let sum = pb.add_row(1.0..=1.0);
    // Which rows each U(h) appears in: U(h) feeds row 2h and row 2 refl^-1(h) + 1.
    let mut inv_refl: HashMap<usize, usize> = HashMap::new();
    for (h, &r) in refl.iter().enumerate() {
        inv_refl.insert(r, h);
    }
    for col in &m.cols {
        let mut factors: Vec<(highs::Row, f64)> = Vec::new();
        for &(r, ups, _) in col {
            if ups == 0 {
                continue;
            }
            let h = r as usize;
            factors.push((rows[2 * h], -(ups as f64)));
            factors.push((rows[2 * inv_refl[&h] + 1], -(ups as f64)));
        }
        factors.push((sum, 2.0));
        pb.add_column(0.0, 0.0.., factors);
    }
    for h in 0..m.rows {
        pb.add_column(1.0, 0.0.., [(rows[2 * h], 1.0), (rows[2 * h + 1], 1.0)]);
    }
    let mut model = pb.optimise(Sense::Maximise);
    configure(&mut model);
    let solved = model.try_solve().ok()?;
    if solved.status() != HighsModelStatus::Optimal {
        return None;
    }
    let sol = solved.get_solution();
    let x = sol.columns()[..m.num_cols()].to_vec();
    let duals = sol.dual_rows();
    let raw: Vec<f64> = (0..m.rows)
        .map(|h| {
            let (pa, pb_) = (duals[2 * h].abs(), duals[2 * h + 1].abs());
            if pa + pb_ > 1e-12 { pa / (pa + pb_) } else { 0.5 }
        })
        .collect();
    // The dual of the symmetric LP only bounds the average over each reflection orbit.
    // Averaging the Labeler with its reflection (write '+' at h as often as '-' at
    // reflect h) makes it a valid Labeler against every single tuple.
    let a = (0..m.rows).map(|h| 0.5 * (raw[h] + 1.0 - raw[refl[h]])).collect();
    Some(OvSolution { value: 1.0 + solved.objective_value(), x, a })
}
