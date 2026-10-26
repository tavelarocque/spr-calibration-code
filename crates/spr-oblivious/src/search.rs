//! eps* search and certification.
//!
//! 1. Bisection on eps with the float feasibility LP of [`crate::lp`]. The interval
//!    `[a, b]` follows the LP verdicts; separately the best feasible strategy seen is
//!    kept with its float eps(pi).
//! 2. Lower bound: the best strategy is rationalised to a common power-of-two
//!    denominator, made an exact distribution, expanded to both reflection members, and
//!    `eps(pi)` is computed exactly ([`crate::exact::verify`]). That exact value is the
//!    certified lower bound.
//! 3. Upper bound: at the infeasible end `b`, a Farkas vector is found by LP, rationalised
//!    and checked exactly against the full unsymmetrised system
//!    ([`crate::exact::farkas_check`]) at the exact rational `b`. If eps* reaches 1/2,
//!    the upper bound 1/2 is the disjointness bound (two later-cell events at the first
//!    prefix are disjoint), whose depth-1 Farkas vector is checked exactly at `1/2 + 2^-60`.

use crate::exact::{farkas_check, verify};
use crate::lp::{Basis, Feas, SymModel, farkas, feasible};
use crate::rational::{dyadic_numerator, exact_f64};
use crate::strategy::Strategy;
use crate::tuples::Space;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
use std::time::Instant;

pub struct Upper {
    pub value: BigRational,
    pub certified: bool,
    pub how: String,
}

pub struct EpsResult {
    pub space: Space,
    pub cols: usize,
    pub rows: usize,
    pub nnz: usize,
    pub lp_solves: usize,
    /// Float bracket from the LP verdicts.
    pub a: f64,
    pub b: f64,
    /// Best feasible symmetric strategy (orbit masses) and its float eps.
    pub x: Vec<f64>,
    pub x_eps: f64,
    pub x_at: f64,
    pub basis: Basis,
    /// Certified lower bound and the exact strategy achieving it.
    pub lower: BigRational,
    pub lower_strategy: Strategy,
    /// The LP basis whose exact basic solution gave the certified lower bound (if the
    /// exact basis polish won).
    pub lower_basis: Option<Basis>,
    pub upper: Upper,
    pub seconds: f64,
}

/// Turn symmetric orbit masses into an exact strategy over all tuples.
pub fn rationalise(m: &SymModel, x: &[f64], bits: u32) -> Strategy {
    let max = x.iter().cloned().fold(0.0, f64::max);
    let mut num: Vec<BigInt> =
        x.iter().map(|&v| if v > 1e-13 * max.max(1e-300) { dyadic_numerator(v, bits) } else { BigInt::zero() }).collect();
    for v in num.iter_mut() {
        if v.is_negative() {
            *v = BigInt::zero();
        }
    }
    // Each orbit has two members, so the masses sum to 2 * sum(num) / 2^bits.
    let target = BigInt::one() << bits;
    let total: BigInt = num.iter().sum::<BigInt>() * 2;
    let diff = &target - &total; // even
    let imax = (0..num.len()).max_by(|&i, &j| num[i].cmp(&num[j])).unwrap();
    num[imax] += diff / 2;
    assert!(!num[imax].is_negative());
    let den = BigInt::one() << bits;
    let mut entries = Vec::new();
    for (o, v) in num.iter().enumerate() {
        if v.is_zero() {
            continue;
        }
        let p = BigRational::new(v.clone(), den.clone());
        let t = m.rep(o).to_vec();
        let rt = m.space.reflect(&t);
        entries.push((t, p.clone()));
        entries.push((rt, p));
    }
    Strategy { n: m.space.n, s: m.space.s, entries }
}

/// Continued-fraction rounding of the orbit masses with denominators up to `max_den`,
/// with the largest mass adjusted so the total is exactly 1.
pub fn rationalise_cf(m: &SymModel, x: &[f64], max_den: u64) -> Option<Strategy> {
    use crate::rational::best_rational;
    let two = BigRational::from_integer(2.into());
    let mut q: Vec<BigRational> =
        x.iter().map(|&v| if v > 1e-12 { best_rational(v, max_den) } else { BigRational::zero() }).collect();
    let total: BigRational = q.iter().cloned().sum::<BigRational>() * &two;
    let imax = (0..q.len()).max_by(|&i, &j| q[i].cmp(&q[j]))?;
    q[imax] += (BigRational::one() - total) / &two;
    if q[imax].is_negative() {
        return None;
    }
    let mut entries = Vec::new();
    for (o, p) in q.iter().enumerate() {
        if p.is_zero() {
            continue;
        }
        let t = m.rep(o).to_vec();
        let rt = m.space.reflect(&t);
        entries.push((t, p.clone()));
        entries.push((rt, p.clone()));
    }
    Some(Strategy { n: m.space.n, s: m.space.s, entries })
}

/// Solve the square basis system exactly at `eps`: nonbasic (tight) rows restricted to
/// basic structural columns, right-hand side 0 on prefix rows and 1 on the sum row.
/// Returns the strategy if every basic value is non-negative.
pub fn basis_solution(m: &SymModel, basis: &Basis, eps: &BigRational) -> Option<Strategy> {
    let cols: Vec<usize> = (0..m.num_cols()).filter(|&o| basis.col_status.get(o) == Some(&1)).collect();
    let rows: Vec<usize> = (0..=m.rows).filter(|&h| basis.row_status.get(h).is_some_and(|&st| st != 1)).collect();
    let k = cols.len();
    if k == 0 || rows.len() != k || k > 400 {
        return None;
    }
    let row_pos: std::collections::HashMap<usize, usize> = rows.iter().enumerate().map(|(i, &r)| (r, i)).collect();
    let mut a = vec![vec![BigRational::zero(); k + 1]; k];
    for (j, &o) in cols.iter().enumerate() {
        for &(row, ups, mem) in &m.cols[o] {
            if let Some(&i) = row_pos.get(&(row as usize)) {
                a[i][j] = BigRational::from_integer(ups.into()) - eps * BigRational::from_integer(mem.into());
            }
        }
        if let Some(&i) = row_pos.get(&m.rows) {
            a[i][j] = BigRational::from_integer(2.into());
        }
    }
    if let Some(&i) = row_pos.get(&m.rows) {
        a[i][k] = BigRational::one();
    }
    // Gaussian elimination with exact pivots.
    for c in 0..k {
        let p = (c..k).find(|&r| !a[r][c].is_zero())?;
        a.swap(c, p);
        let piv = a[c][c].clone();
        for v in a[c].iter_mut() {
            *v = &*v / &piv;
        }
        let pivot_row = a[c].clone();
        for (r, row) in a.iter_mut().enumerate() {
            if r != c && !row[c].is_zero() {
                let f = row[c].clone();
                for (v, p) in row.iter_mut().zip(pivot_row.iter()).skip(c) {
                    *v -= &f * p;
                }
            }
        }
    }
    let mut entries = Vec::new();
    for (j, &o) in cols.iter().enumerate() {
        let v = a[j][k].clone();
        if v.is_negative() {
            return None;
        }
        if v.is_zero() {
            continue;
        }
        let t = m.rep(o).to_vec();
        let rt = m.space.reflect(&t);
        entries.push((t, v.clone()));
        entries.push((rt, v));
    }
    Some(Strategy { n: m.space.n, s: m.space.s, entries })
}

/// Largest rational eps' (by exact bisection between `at - 1e-6` and `b + 1e-9`) at
/// which the basis gives a strategy with exact eps(pi) >= eps'. Returns that strategy
/// and its exact eps.
fn polish(m: &SymModel, basis: &Basis, at: f64, b: f64) -> Option<(Strategy, BigRational)> {
    let ok = |e: &BigRational| -> Option<(Strategy, BigRational)> {
        let st = basis_solution(m, basis, e)?;
        let v = verify(&st).ok()?.eps;
        (v >= *e).then_some((st, v))
    };
    let mut lo = exact_f64((at - 1e-6).max(0.0));
    let mut best = ok(&lo)?;
    let mut hi = exact_f64(b + 1e-9);
    if let Some(r) = ok(&hi) {
        return Some(r);
    }
    let two = BigRational::from_integer(2.into());
    for _ in 0..50 {
        let mid = (&lo + &hi) / &two;
        match ok(&mid) {
            Some(r) => {
                lo = mid;
                if r.1 > best.1 {
                    best = r;
                }
            }
            None => hi = mid,
        }
    }
    Some(best)
}

/// Exact Farkas vector at `eps` from the float Farkas LP, or None.
fn farkas_at(m: &SymModel, eps: &BigRational) -> bool {
    let e = to_f64_r(eps);
    let Some(y) = farkas(m, e) else { return false };
    let ymax = y.iter().cloned().fold(0.0, f64::max).max(1e-300);
    let yr: Vec<BigRational> = y
        .iter()
        .map(|&v| {
            let v = if v > 1e-14 * ymax { v / ymax } else { 0.0 };
            BigRational::new(dyadic_numerator(v, 50), BigInt::one() << 50u32)
        })
        .collect();
    farkas_check(&m.space, &yr, eps).is_ok()
}

fn to_f64_r(r: &BigRational) -> f64 {
    crate::rational::to_f64(r)
}

/// Push a certified upper bound down towards the certified lower bound by exact
/// bisection: each candidate needs a Farkas vector that passes the exact check.
fn tighten_upper(m: &SymModel, lower: &BigRational, upper: &mut Upper) {
    let two = BigRational::from_integer(2.into());
    let mut lo = lower.clone();
    let mut hi = upper.value.clone();
    for _ in 0..25 {
        if to_f64_r(&(&hi - &lo)) < 1e-12 {
            break;
        }
        let mid = (&lo + &hi) / &two;
        // Round mid to a short dyadic so the check stays cheap.
        let mid = exact_f64(to_f64_r(&mid));
        if mid <= lo || mid >= hi {
            break;
        }
        if farkas_at(m, &mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    if hi < upper.value {
        upper.value = hi;
        upper.how = "rational Farkas vector checked exactly (tightened by exact bisection)".into();
    }
}

fn certify_upper(m: &SymModel, b: f64) -> Upper {
    if b >= 0.5 {
        // Depth-1 rows with weight 1: every column has value up + down - 2 eps <= 1 - 2 eps.
        let sp = m.space;
        let mut y = vec![BigRational::zero(); sp.num_prefixes()];
        for v in y.iter_mut().take(sp.n) {
            *v = BigRational::one();
        }
        let eps = BigRational::new(BigInt::one(), BigInt::from(2)) + BigRational::new(BigInt::one(), BigInt::one() << 60);
        return match farkas_check(&sp, &y, &eps) {
            Ok(_) => Upper {
                value: BigRational::new(1.into(), 2.into()),
                certified: true,
                how: "disjointness bound; depth-1 Farkas vector checked exactly at 1/2 + 2^-60".into(),
            },
            Err(e) => Upper { value: exact_f64(0.5), certified: false, how: format!("trivial certificate failed: {e}") },
        };
    }
    let mut eps = b;
    for attempt in 0..4 {
        if let Some(y) = farkas(m, eps) {
            let ymax = y.iter().cloned().fold(0.0, f64::max).max(1e-300);
            let yr: Vec<BigRational> = y
                .iter()
                .map(|&v| {
                    let v = if v > 1e-14 * ymax { v / ymax } else { 0.0 };
                    BigRational::new(dyadic_numerator(v, 50), BigInt::one() << 50u32)
                })
                .collect();
            let e = exact_f64(eps);
            if farkas_check(&m.space, &yr, &e).is_ok() {
                return Upper {
                    value: e,
                    certified: true,
                    how: format!("rational Farkas vector checked exactly (attempt {attempt})"),
                };
            }
        }
        eps += 1e-9 * (attempt + 1) as f64;
    }
    Upper { value: exact_f64(b), certified: false, how: "float, uncertified".into() }
}

/// Largest number of tuples the LP is built for (about 3 GB at the measured 200 bytes per
/// nonzero); larger instances are refused rather than left to run out of memory.
pub const MAX_TUPLES: usize = 5_000_000;

pub fn solve(space: Space, tol: f64) -> Result<EpsResult, String> {
    if space.num_tuples() > MAX_TUPLES {
        return Err(format!("{} tuples exceeds the limit of {MAX_TUPLES}", space.num_tuples()));
    }
    let t0 = Instant::now();
    let m = SymModel::build(space);
    let mut solves = 0;
    let mut best: Option<(Vec<f64>, f64, f64, Basis)> = None; // (x, eps(x), eps used, basis)
    let record = |x: Vec<f64>, at: f64, basis: Basis, best: &mut Option<(Vec<f64>, f64, f64, Basis)>| {
        let e = m.eps_of(&x);
        if best.as_ref().is_none_or(|b| e > b.1) {
            *best = Some((x, e, at, basis));
        }
    };
    let (mut a, mut b);
    solves += 1;
    match feasible(&m, 0.5) {
        Feas::Feasible { x, basis } => {
            record(x, 0.5, basis, &mut best);
            a = 0.5;
            b = 0.5;
        }
        Feas::Infeasible => {
            a = 0.0;
            b = 0.5;
        }
        Feas::Unknown(s) => return Err(format!("HiGHS status {s} at eps = 1/2")),
    }
    if best.is_none() {
        solves += 1;
        match feasible(&m, 0.0) {
            Feas::Feasible { x, basis } => record(x, 0.0, basis, &mut best),
            _ => return Err("LP infeasible at eps = 0".into()),
        }
    }
    let mut unknown_stop = false;
    let mut bases: Vec<(f64, Basis)> = Vec::new();
    while b - a > tol {
        let mid = 0.5 * (a + b);
        solves += 1;
        match feasible(&m, mid) {
            Feas::Feasible { x, basis } => {
                bases.push((mid, basis.clone()));
                record(x, mid, basis, &mut best);
                a = mid;
            }
            Feas::Infeasible => b = mid,
            Feas::Unknown(_) => {
                // Undecided even after the retries: an exactly verified Farkas vector at
                // mid settles it; otherwise stop with the bracket known so far.
                if certify_upper(&m, mid).certified {
                    b = mid;
                } else {
                    unknown_stop = true;
                    break;
                }
            }
        }
    }
    let (x, x_eps, x_at, basis) = best.unwrap();
    // Candidates: dyadic rounding, and continued-fraction rounding at several
    // denominator bounds (which recovers small-denominator strategies exactly). The
    // exactly verified best one is the certified lower bound.
    // Masses below a threshold are LP noise that creates prefixes with tiny mu, where
    // the ratio U/mu is badly conditioned; candidates with those masses removed are
    // tried as well.
    let mut lower_strategy = rationalise(&m, &x, 50);
    let mut lower = verify(&lower_strategy)?.eps;
    for cut in [0.0, 1e-9, 1e-7, 1e-5] {
        let xc: Vec<f64> = x.iter().map(|&v| if v < cut { 0.0 } else { v }).collect();
        let mut candidates = vec![rationalise(&m, &xc, 50)];
        for max_den in [100u64, 10_000, 1_000_000, 100_000_000] {
            candidates.extend(rationalise_cf(&m, &xc, max_den));
        }
        for st in candidates {
            let e = verify(&st)?.eps;
            if e > lower {
                lower = e;
                lower_strategy = st;
            }
        }
    }
    // Exact basis polish: for the bases of the last feasible solves, the basic solution
    // is computed exactly at rational eps', and eps' is pushed up by exact bisection
    // while that solution stays a valid strategy with eps(pi) >= eps'.
    bases.push((x_at, basis.clone()));
    bases.sort_by(|p, q| q.0.partial_cmp(&p.0).unwrap());
    bases.dedup_by(|p, q| p.1.col_status == q.1.col_status && p.1.row_status == q.1.row_status);
    let mut lower_basis = None;
    for (at, bs) in bases.iter().take(6) {
        if let Some((st, e)) = polish(&m, bs, *at, b)
            && e > lower
        {
            lower = e;
            lower_strategy = st;
            lower_basis = Some(bs.clone());
        }
    }
    let mut upper = certify_upper(&m, b);
    if upper.certified && b < 0.5 {
        tighten_upper(&m, &lower, &mut upper);
    }
    if unknown_stop {
        eprintln!("note: bisection stopped at a HiGHS Unknown status; bracket [{a}, {b}]");
    }
    Ok(EpsResult {
        space,
        cols: m.num_cols(),
        rows: m.rows,
        nnz: m.nnz(),
        lp_solves: solves,
        a,
        b,
        x,
        x_eps,
        x_at,
        basis,
        lower,
        lower_strategy,
        lower_basis,
        upper,
        seconds: t0.elapsed().as_secs_f64(),
    })
}

pub struct OvResult {
    pub float: f64,
    /// Certified bounds: lower from an exact strategy, upper from an exact Labeler.
    pub lower: BigRational,
    pub upper: BigRational,
    pub strategy: Strategy,
    pub seconds: f64,
}

/// ov* with exact bounds. The LP primal (strategy) and dual (randomised Labeler) are
/// rationalised by continued fractions at several denominator bounds; each candidate is
/// evaluated exactly, and the best lower and upper bounds are kept. ov* is certified
/// exactly when they coincide.
pub fn ov_certified(space: Space) -> Result<OvResult, String> {
    use crate::exact::ov_upper;
    if space.num_tuples() > MAX_TUPLES {
        return Err(format!("{} tuples exceeds the limit of {MAX_TUPLES}", space.num_tuples()));
    }
    use crate::rational::best_rational;
    let t0 = Instant::now();
    let m = SymModel::build(space);
    let sol = crate::lp::ov(&m).ok_or("ov LP not solved to optimality")?;
    let mut best_lo: Option<(BigRational, Strategy)> = None;
    let mut best_hi: Option<BigRational> = None;
    for max_den in [100u64, 10_000, 1_000_000, 100_000_000, 10_000_000_000] {
        // Primal: orbit masses with bounded denominators, then make the total exactly 1.
        let mut q: Vec<BigRational> =
            sol.x.iter().map(|&v| if v > 1e-12 { best_rational(v, max_den) } else { BigRational::zero() }).collect();
        let total: BigRational = q.iter().cloned().sum::<BigRational>() * BigRational::from_integer(2.into());
        let imax = (0..q.len()).max_by(|&i, &j| q[i].cmp(&q[j])).unwrap();
        q[imax] += (BigRational::one() - total) / BigRational::from_integer(2.into());
        if !q[imax].is_negative() {
            let mut entries = Vec::new();
            for (o, p) in q.iter().enumerate() {
                if p.is_zero() {
                    continue;
                }
                let t = m.rep(o).to_vec();
                let rt = space.reflect(&t);
                entries.push((t, p.clone()));
                entries.push((rt, p.clone()));
            }
            let st = Strategy { n: space.n, s: space.s, entries };
            let ov = verify(&st)?.ov;
            if best_lo.as_ref().is_none_or(|b| ov > b.0) {
                best_lo = Some((ov, st));
            }
        }
        // Dual: Labeler probabilities in [0, 1].
        let a: Vec<BigRational> = sol
            .a
            .iter()
            .map(|&v| {
                let r = best_rational(v.clamp(0.0, 1.0), max_den);
                r.max(BigRational::zero()).min(BigRational::one())
            })
            .collect();
        let up = ov_upper(&space, &a)?;
        if best_hi.as_ref().is_none_or(|b| &up < b) {
            best_hi = Some(up);
        }
    }
    let (lower, strategy) = best_lo.ok_or("no rational primal found")?;
    Ok(OvResult { float: sol.value, lower, upper: best_hi.unwrap(), strategy, seconds: t0.elapsed().as_secs_f64() })
}
