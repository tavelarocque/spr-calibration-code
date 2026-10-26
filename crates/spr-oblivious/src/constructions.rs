//! Explicit strategies: the DDFGKO binomial cube, oblivious binary search, and the
//! product of two strategies.

use crate::strategy::Strategy;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use std::collections::HashMap;

/// All strings over `alphabet` of length `d` with exactly `zeros` zero entries, in
/// lexicographic order (the alphabet is given in increasing order).
fn strings_with_zeros(alphabet: &[i8], d: usize, zeros: usize) -> Vec<Vec<i8>> {
    let mut out = Vec::new();
    let mut cur = Vec::with_capacity(d);
    fn rec(alpha: &[i8], d: usize, zeros: usize, cur: &mut Vec<i8>, out: &mut Vec<Vec<i8>>) {
        let z = cur.iter().filter(|&&c| c == 0).count();
        if cur.len() == d {
            if z == zeros {
                out.push(cur.clone());
            }
            return;
        }
        for &a in alpha {
            if a == 0 && z == zeros {
                continue;
            }
            cur.push(a);
            rec(alpha, d, zeros, cur, out);
            cur.pop();
        }
    }
    rec(alphabet, d, zeros, &mut cur, &mut out);
    out
}

/// The binomial cube of Dagan et al. (arXiv:2406.13668), Lemma 7.1 in version 4
/// (Lemma 6.1 in version 2).
///
/// Time labels: binary strings `w` of length `d` with exactly `k` zeros, in
/// lexicographic order. Cells: strings `q` in `{-1, 0, 1}^d` with exactly `k` zeros,
/// ordered lexicographically with `-1 < 0 < 1`, numbered 0, 1, .... One fair coin
/// `xi_u` per binary prefix `u`; step `t` points at `q^(t)_l = w^(t)_l * xi_{w^(t)_1..l-1}`.
/// Every outcome of the coins that are actually used gets equal mass.
pub fn binomial_cube(d: usize, k: usize) -> Result<Strategy, String> {
    assert!(k <= d);
    let labels: Vec<Vec<i8>> = strings_with_zeros(&[0, 1], d, k);
    let cells: Vec<Vec<i8>> = strings_with_zeros(&[-1, 0, 1], d, k);
    if cells.len() > 255 {
        return Err(format!("{} cells do not fit in u8", cells.len()));
    }
    let cell_index: HashMap<Vec<i8>, u8> = cells.iter().enumerate().map(|(i, q)| (q.clone(), i as u8)).collect();
    // Coins that are ever multiplied by a 1: prefixes w[..l] with w[l] = 1.
    let mut coins: Vec<Vec<i8>> = Vec::new();
    for w in &labels {
        for l in 0..d {
            if w[l] == 1 {
                coins.push(w[..l].to_vec());
            }
        }
    }
    coins.sort();
    coins.dedup();
    let c = coins.len();
    if c > 24 {
        return Err(format!("{c} coins: 2^{c} outcomes is too many to enumerate"));
    }
    let coin_index: HashMap<Vec<i8>, usize> = coins.iter().enumerate().map(|(i, u)| (u.clone(), i)).collect();
    let mass = BigRational::new(BigInt::one(), BigInt::one() << c);
    let mut entries = Vec::with_capacity(1 << c);
    for outcome in 0u64..(1u64 << c) {
        let xi = |u: &[i8]| -> i8 { if outcome >> coin_index[u] & 1 == 1 { 1 } else { -1 } };
        let tuple: Vec<u8> = labels
            .iter()
            .map(|w| {
                let q: Vec<i8> = (0..d).map(|l| if w[l] == 1 { xi(&w[..l]) } else { 0 }).collect();
                cell_index[&q]
            })
            .collect();
        entries.push((tuple, mass.clone()));
    }
    Ok(Strategy { n: cells.len(), s: labels.len(), entries })
}

/// Number of coins the binomial cube uses (for reporting).
pub fn binomial_cube_coins(d: usize, k: usize) -> usize {
    let labels = strings_with_zeros(&[0, 1], d, k);
    let mut coins: Vec<Vec<i8>> = Vec::new();
    for w in &labels {
        for l in 0..d {
            if w[l] == 1 {
                coins.push(w[..l].to_vec());
            }
        }
    }
    coins.sort();
    coins.dedup();
    coins.len()
}

/// Oblivious binary search on the first `2^s - 1` of `n` cells: point at the middle of
/// the current interval, then continue in a uniformly random half.
pub fn binary_search(n: usize, s: usize) -> Strategy {
    assert!(n >= (1 << s) - 1);
    let mut entries = Vec::new();
    for path in 0u64..(1 << (s - 1)) {
        let (mut lo, mut hi) = (0i64, (1i64 << s) - 2);
        let mut t = Vec::with_capacity(s);
        for step in 0..s {
            let mid = (lo + hi) / 2;
            t.push(mid as u8);
            if step + 1 < s {
                if path >> step & 1 == 1 { lo = mid + 1 } else { hi = mid - 1 }
            }
        }
        entries.push((t, BigRational::new(BigInt::one(), BigInt::one() << (s - 1))));
    }
    Strategy { n, s, entries }
}

/// The product of two strategies. Cells are pairs `(a, b)`
/// ordered lexicographically, numbered `a * n2 + b`. Draw an outer tuple from `p1`;
/// for each outer cell in turn draw a fresh independent inner tuple from `p2` and
/// point at `(a, b_1), ..., (a, b_{s2})`.
pub fn product(p1: &Strategy, p2: &Strategy) -> Result<Strategy, String> {
    let n = p1.n * p2.n;
    if n > 255 {
        return Err(format!("product has {n} cells, more than u8 allows"));
    }
    let s = p1.s * p2.s;
    let mut entries = Vec::new();
    for (outer, po) in &p1.entries {
        // Mixed-radix enumeration of one inner tuple per outer cell.
        let m = p2.entries.len();
        let combos = m.checked_pow(p1.s as u32).ok_or("too many combinations")?;
        for mut code in 0..combos {
            let mut tuple = Vec::with_capacity(s);
            let mut mass = po.clone();
            for &a in outer {
                let (inner, pi) = &p2.entries[code % m];
                code /= m;
                mass *= pi;
                tuple.extend(inner.iter().map(|&b| (a as usize * p2.n + b as usize) as u8));
            }
            entries.push((tuple, mass));
        }
    }
    Ok(Strategy { n, s, entries })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::exact::{verify, verify_naive};

    #[test]
    fn small_cube_by_hand() {
        // d = 2, k = 1: cells (-1,0) < (0,-1) < (0,1) < (1,0); step 1 is cell 1 or 2,
        // step 2 is cell 0 or 3, independently.
        let st = binomial_cube(2, 1).unwrap();
        assert_eq!((st.n, st.s), (4, 2));
        let mut tuples: Vec<Vec<u8>> = st.entries.iter().map(|e| e.0.clone()).collect();
        tuples.sort();
        assert_eq!(tuples, vec![vec![1, 0], vec![1, 3], vec![2, 0], vec![2, 3]]);
        assert_eq!(verify(&st).unwrap().eps, BigRational::new(1.into(), 2.into()));
    }

    #[test]
    fn binary_search_is_one_half() {
        for s in 2..=4 {
            let st = binary_search((1 << s) - 1, s);
            let half = BigRational::new(1.into(), 2.into());
            assert_eq!(verify(&st).unwrap().eps, half);
            assert_eq!(verify_naive(&st).unwrap(), half);
        }
    }
}
