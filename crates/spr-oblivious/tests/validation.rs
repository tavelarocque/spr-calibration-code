//! Validation tests for the oblivious solver.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};
use spr_oblivious::constructions::{binary_search, binomial_cube, product};
use spr_oblivious::exact::{verify, verify_naive};
use spr_oblivious::search::solve;
use spr_oblivious::strategy::Strategy;
use spr_oblivious::tuples::Space;

fn r(a: i64, b: i64) -> BigRational {
    BigRational::new(a.into(), b.into())
}

fn half() -> BigRational {
    r(1, 2)
}

/// (a) eps*(n, 2) = 1/2 for n >= 3, certified on both sides.
#[test]
fn a_eps_n_2_is_one_half() {
    for n in 3..=12 {
        let res = solve(Space::new(n, 2), 1e-10).unwrap();
        assert_eq!(res.lower, half(), "n={n}");
        assert!(res.upper.certified && res.upper.value == half(), "n={n}");
    }
    // n = 2: the two cells cannot give 1/2 (the first sign is always erasable on one side).
    let res = solve(Space::new(2, 2), 1e-10).unwrap();
    assert!(res.lower < half());
}

/// (a) eps*(15, 4) >= 1/2 from the oblivious binary search, verified exactly.
#[test]
fn a_binary_search_15_4() {
    let st = binary_search(15, 4);
    assert_eq!(verify(&st).unwrap().eps, half());
    assert_eq!(verify_naive(&st).unwrap(), half());
}

/// (a) eps <= 1/2 for every s >= 2: every certified upper bound on a small grid is at
/// most 1/2, and every certified lower bound is too.
#[test]
fn a_never_above_one_half() {
    for s in 2..=4 {
        for n in s..=9 {
            let res = solve(Space::new(n, s), 1e-9).unwrap();
            assert!(res.upper.certified, "n={n} s={s}");
            assert!(res.upper.value <= half() + r(1, 1 << 30) && res.lower <= half(), "n={n} s={s}");
        }
    }
}

/// (d) The product lemma, verified exactly on (3,2) x (3,2) and (3,2) x (7,3).
#[test]
fn d_product_lemma() {
    let p32 = solve(Space::new(3, 2), 1e-10).unwrap().lower_strategy;
    let p73 = binary_search(7, 3);
    for (a, b) in [(&p32, &p32), (&p32, &p73), (&p73, &p32)] {
        let (ea, eb) = (verify(a).unwrap().eps, verify(b).unwrap().eps);
        let prod = product(a, b).unwrap();
        let e = verify(&prod).unwrap().eps;
        assert!(e >= &ea * &eb, "product of ({},{}) and ({},{}): {e} < {}", a.n, a.s, b.n, b.s, &ea * &eb);
        assert_eq!(verify_naive(&prod).unwrap(), e);
    }
}

fn perturbed(st: &Strategy, from: usize, to: usize, delta: &BigRational) -> Strategy {
    let mut s = st.clone();
    s.entries[from].1 -= delta;
    s.entries[to].1 += delta;
    s
}

/// (e) The verifier rejects perturbed strategies exactly when the perturbation drops eps
/// below the claim. Every verdict is cross-checked against the independent verifier.
#[test]
fn e_verifier_rejects_broken_strategies() {
    // A deliberate failure: (3, 2) with the end cells at 0.6 / 0.4.
    let st = Strategy { n: 3, s: 2, entries: vec![(vec![1, 0], r(3, 5)), (vec![1, 2], r(2, 5))] };
    let e = verify(&st).unwrap().eps;
    assert_eq!(e, r(2, 5));
    assert!(e < half());
    // Masses that do not sum to 1, a repeated cell, an out-of-range cell: rejected outright.
    assert!(verify(&Strategy { n: 3, s: 2, entries: vec![(vec![1, 0], r(1, 2))] }).is_err());
    assert!(verify(&Strategy { n: 3, s: 2, entries: vec![(vec![1, 1], r(1, 1))] }).is_err());
    assert!(verify(&Strategy { n: 3, s: 2, entries: vec![(vec![1, 3], r(1, 1))] }).is_err());

    // Systematic: every single mass transfer of size 1/1000 between support tuples of the
    // certified (6, 3) strategy.
    let base = solve(Space::new(6, 3), 1e-10).unwrap();
    let claim = base.lower.clone();
    let st = base.lower_strategy;
    let delta = r(1, 1000);
    let (mut rejected, mut kept) = (0, 0);
    for i in 0..st.entries.len() {
        for j in 0..st.entries.len() {
            if i == j || st.entries[i].1 < delta {
                continue;
            }
            let p = perturbed(&st, i, j, &delta);
            let fast = verify(&p).unwrap().eps;
            let slow = verify_naive(&p).unwrap();
            assert_eq!(fast, slow);
            if fast < claim {
                rejected += 1;
            } else {
                kept += 1;
            }
        }
    }
    assert!(rejected > 0, "no perturbation lowered eps");
    eprintln!("(6,3) perturbations: {rejected} rejected, {kept} kept the claim");
}

/// Addition A: the DDFGKO binomial cube has eps = 2^-k exactly.
#[test]
fn addition_a_binomial_cube() {
    for (d, k) in [(3, 1), (4, 1), (4, 2), (5, 2)] {
        let st = binomial_cube(d, k).unwrap();
        let want = BigRational::new(BigInt::one(), BigInt::one() << k);
        assert_eq!(verify(&st).unwrap().eps, want, "d={d} k={k}");
        if st.entries.len() <= 5000 {
            assert_eq!(verify_naive(&st).unwrap(), want, "d={d} k={k}");
        }
    }
}

/// eps* is non-decreasing in n and non-increasing in s (each by a short argument). Checked on certified intervals:
/// lower(n, s) <= upper(n + 1, s) and lower(n, s + 1) <= upper(n, s).
#[test]
fn monotonicity_in_n_and_s() {
    let mut lo = std::collections::HashMap::new();
    let mut hi = std::collections::HashMap::new();
    for s in 2..=5 {
        for n in s..=9 {
            let res = solve(Space::new(n, s), 1e-10).unwrap();
            assert!(res.upper.certified);
            lo.insert((n, s), res.lower);
            hi.insert((n, s), res.upper.value);
        }
    }
    for (&(n, s), l) in &lo {
        if let Some(u) = hi.get(&(n + 1, s)) {
            assert!(l <= u, "n monotonicity fails at ({n},{s})");
        }
        if let Some(l2) = lo.get(&(n, s + 1)) {
            assert!(l2 <= &hi[&(n, s)], "s monotonicity fails at ({n},{s})");
        }
    }
    // eps* = 0 is genuine when n is small for s (for example n = s), so only
    // non-negativity is required.
    assert!(!lo.is_empty() && lo.values().all(|v| v >= &BigRational::zero()));
}
