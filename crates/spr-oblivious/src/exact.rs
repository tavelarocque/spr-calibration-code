//! Exact certification. Everything here is BigInt / BigRational; nothing in this
//! module touches a float.
//!
//! All checks use the full, unsymmetrised definitions:
//! `mu(h)`, `U(h)` (every later cell strictly greater than the last cell of `h`) and
//! `D(h)` (every later cell strictly smaller), over every prefix `h` of length
//! `1 ..= s - 1`. Prefixes of length `s` have `U = D = mu` and never bind.

use crate::strategy::Strategy;
use crate::tuples::{Space, up_down};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
use std::collections::BTreeMap;

/// The result of an exact epsilon computation.
#[derive(Clone, Debug)]
pub struct EpsReport {
    /// `min over prefixes h with mu(h) > 0 of min(U(h), D(h)) / mu(h)`; 1 if `s = 1`.
    pub eps: BigRational,
    /// A prefix attaining the minimum, and which side ('U' or 'D').
    pub worst: Option<(Vec<u8>, char)>,
    /// Number of prefixes with positive mass.
    pub prefixes: usize,
    /// `ov(pi) = sum over all prefixes of length 1..=s of min(U, D)`.
    pub ov: BigRational,
}

/// Sorted support with integer masses, and their common denominator.
type IntegerMasses = (Vec<(Vec<u8>, BigInt)>, BigInt);

/// Validate a strategy and return it as sorted tuples with integer masses over a
/// common denominator.
fn normalise(st: &Strategy) -> Result<IntegerMasses, String> {
    let mut merged: BTreeMap<Vec<u8>, BigRational> = BTreeMap::new();
    for (t, p) in &st.entries {
        if t.len() != st.s {
            return Err(format!("tuple {t:?} has length {} not {}", t.len(), st.s));
        }
        let mut seen = 0u128;
        for &c in t {
            if c as usize >= st.n {
                return Err(format!("cell {c} out of range in {t:?}"));
            }
            if seen >> c & 1 == 1 {
                return Err(format!("tuple {t:?} repeats cell {c}"));
            }
            seen |= 1u128 << c;
        }
        if p.is_negative() {
            return Err(format!("negative probability for {t:?}"));
        }
        *merged.entry(t.clone()).or_insert_with(BigRational::zero) += p;
    }
    let total: BigRational = merged.values().cloned().sum();
    if total != BigRational::one() {
        return Err(format!("probabilities sum to {total}, not 1"));
    }
    let lcm = merged.values().fold(BigInt::one(), |acc, p| acc.lcm(p.denom()));
    let out: Vec<(Vec<u8>, BigInt)> = merged
        .into_iter()
        .filter(|(_, p)| !p.is_zero())
        .map(|(t, p)| (t, (p * BigRational::from_integer(lcm.clone())).to_integer()))
        .collect();
    Ok((out, lcm))
}

/// Exact `eps(pi)` and `ov(pi)` by one sweep per prefix length over the sorted support.
pub fn verify(st: &Strategy) -> Result<EpsReport, String> {
    let (tuples, lcm) = normalise(st)?;
    let flags: Vec<(Vec<bool>, Vec<bool>)> = tuples.iter().map(|(t, _)| up_down(t)).collect();
    let mut best: Option<(BigInt, BigInt, Vec<u8>, char)> = None; // (num, den) of the min ratio
    let mut prefixes = 0usize;
    let mut ov_num = lcm.clone(); // the length-s prefixes contribute total mass 1
    for i in 1..st.s {
        let mut a = 0;
        while a < tuples.len() {
            let mut b = a;
            let (mut mu, mut u, mut d) = (BigInt::zero(), BigInt::zero(), BigInt::zero());
            while b < tuples.len() && tuples[b].0[..i] == tuples[a].0[..i] {
                mu += &tuples[b].1;
                if flags[b].0[i - 1] {
                    u += &tuples[b].1;
                }
                if flags[b].1[i - 1] {
                    d += &tuples[b].1;
                }
                b += 1;
            }
            if mu.is_positive() {
                prefixes += 1;
                ov_num += u.clone().min(d.clone());
                for (side, val) in [('U', &u), ('D', &d)] {
                    let better = match &best {
                        None => true,
                        Some((bn, bd, _, _)) => val * bd < bn * &mu,
                    };
                    if better {
                        best = Some((val.clone(), mu.clone(), tuples[a].0[..i].to_vec(), side));
                    }
                }
            }
            a = b;
        }
    }
    let ov = BigRational::new(ov_num, lcm);
    Ok(match best {
        None => EpsReport { eps: BigRational::one(), worst: None, prefixes, ov },
        Some((num, den, p, side)) => EpsReport { eps: BigRational::new(num, den), worst: Some((p, side)), prefixes, ov },
    })
}

/// A second, independent exact implementation, written straight from Definition 5.1
/// with no sorting, no sweep and no common denominator: for every prefix of the
/// support, sum conditional masses tuple by tuple. Quadratic; for small supports and
/// for re-verifying any surprising result.
pub fn verify_naive(st: &Strategy) -> Result<BigRational, String> {
    let total: BigRational = st.entries.iter().map(|(_, p)| p.clone()).sum();
    if total != BigRational::one() || st.entries.iter().any(|(_, p)| p.is_negative()) {
        return Err("not a probability distribution".into());
    }
    for (t, _) in &st.entries {
        let mut v = t.clone();
        v.sort();
        v.dedup();
        if v.len() != st.s || t.iter().any(|&c| c as usize >= st.n) {
            return Err(format!("bad tuple {t:?}"));
        }
    }
    let mut prefixes: Vec<Vec<u8>> = Vec::new();
    for (t, p) in &st.entries {
        if p.is_zero() {
            continue;
        }
        for i in 1..st.s {
            prefixes.push(t[..i].to_vec());
        }
    }
    prefixes.sort();
    prefixes.dedup();
    let mut eps = BigRational::one();
    for h in &prefixes {
        let i = h.len();
        let last = h[i - 1];
        let mut mu = BigRational::zero();
        let mut up = BigRational::zero();
        let mut down = BigRational::zero();
        for (t, p) in &st.entries {
            if t[..i] != h[..] {
                continue;
            }
            mu += p;
            if t[i..].iter().all(|&c| c > last) {
                up += p;
            }
            if t[i..].iter().all(|&c| c < last) {
                down += p;
            }
        }
        if mu.is_positive() {
            eps = eps.min(&up / &mu).min(&down / &mu);
        }
    }
    Ok(eps)
}

/// Farkas certificate check for "no pi has eps(pi) >= eps".
///
/// The primal system is: `U(h) - eps mu(h) >= 0` and `D(h) - eps mu(h) >= 0` for every
/// prefix `h`, `sum pi = 1`, `pi >= 0`. A vector `y >= 0` on these rows with
/// `y^T A < 0` in every column (every tuple) proves infeasibility: for a feasible pi,
/// `y^T A pi >= 0`, but `y^T A pi` is a convex combination of negative numbers.
///
/// `y_u[idx]` weights the U-row of the prefix with dense index `idx`; the D-row of `h`
/// gets `y_u` of the reflected prefix (the expansion of a certificate found for the
/// symmetric LP). The check runs over every tuple of `[n]^s` with the full U and D rows.
/// Returns the largest column value, which must be negative.
pub fn farkas_check(space: &Space, y_u: &[BigRational], eps: &BigRational) -> Result<BigRational, String> {
    if y_u.len() != space.num_prefixes() {
        return Err("certificate has the wrong length".into());
    }
    if y_u.iter().any(|y| y.is_negative()) {
        return Err("certificate has a negative entry".into());
    }
    let lcm = y_u.iter().fold(BigInt::one(), |acc, y| acc.lcm(y.denom()));
    let y: Vec<BigInt> = y_u.iter().map(|v| (v * BigRational::from_integer(lcm.clone())).to_integer()).collect();
    let (p, q) = (eps.numer().clone(), eps.denom().clone());
    let all = space.all_tuples();
    let mut worst: Option<BigInt> = None;
    for t in all.chunks(space.s) {
        let (up, down) = up_down(t);
        // column value times lcm * q
        let mut val = BigInt::zero();
        for i in 1..space.s {
            let h = &t[..i];
            let yu = &y[space.prefix_index(h)];
            let yd = &y[space.prefix_index(&space.reflect(h))];
            let cu = if up[i - 1] { &q - &p } else { -&p };
            let cd = if down[i - 1] { &q - &p } else { -&p };
            val += yu * cu + yd * cd;
        }
        if !val.is_negative() {
            return Err(format!("column {t:?} is not negative: {val}"));
        }
        if worst.as_ref().is_none_or(|w| &val > w) {
            worst = Some(val);
        }
    }
    Ok(BigRational::new(worst.unwrap_or_else(BigInt::zero), lcm * q))
}

/// Exact upper bound on ov* from a randomised Labeler: at prefix `h` the Labeler writes
/// '+' with probability `a[idx(h)]` (and '-' otherwise). For every pi,
/// `ov(pi) <= 1 + sum_h (a_h U(h) + (1 - a_h) D(h))`, whose maximum over pi is attained
/// at a single tuple. Returns that maximum over all tuples of `[n]^s`.
pub fn ov_upper(space: &Space, a: &[BigRational]) -> Result<BigRational, String> {
    if a.len() != space.num_prefixes() {
        return Err("Labeler strategy has the wrong length".into());
    }
    if a.iter().any(|x| x.is_negative() || x > &BigRational::one()) {
        return Err("Labeler probabilities must lie in [0, 1]".into());
    }
    let lcm = a.iter().fold(BigInt::one(), |acc, y| acc.lcm(y.denom()));
    let an: Vec<BigInt> = a.iter().map(|v| (v * BigRational::from_integer(lcm.clone())).to_integer()).collect();
    let all = space.all_tuples();
    let mut best = BigInt::zero();
    for t in all.chunks(space.s) {
        let (up, down) = up_down(t);
        let mut val = lcm.clone();
        for i in 1..space.s {
            let ah = &an[space.prefix_index(&t[..i])];
            if up[i - 1] {
                val += ah;
            }
            if down[i - 1] {
                val += &lcm - ah;
            }
        }
        if val > best {
            best = val;
        }
    }
    Ok(BigRational::new(best, lcm))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(a: i64, b: i64) -> BigRational {
        BigRational::new(a.into(), b.into())
    }

    #[test]
    fn three_two_strategy() {
        // Point at the middle cell, then at a uniformly random end.
        let st = Strategy { n: 3, s: 2, entries: vec![(vec![1, 0], r(1, 2)), (vec![1, 2], r(1, 2))] };
        let rep = verify(&st).unwrap();
        assert_eq!(rep.eps, r(1, 2));
        assert_eq!(verify_naive(&st).unwrap(), r(1, 2));
        // ov = 1 (last sign) + min(1/2, 1/2)
        assert_eq!(rep.ov, r(3, 2));
    }

    #[test]
    fn rejects_bad_input() {
        let st = Strategy { n: 3, s: 2, entries: vec![(vec![1, 1], r(1, 1))] };
        assert!(verify(&st).is_err());
        let st = Strategy { n: 3, s: 2, entries: vec![(vec![1, 0], r(1, 2))] };
        assert!(verify(&st).is_err());
    }
}
