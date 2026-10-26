//! Validation (b) on random reachable states, (c) proptest invariants, and (d) the
//! Qiao-Valiant known value.

use proptest::prelude::*;
use spr_core::{Board, Forced, Full, Qv, Variant};
use spr_solve::brute::Brute;
use spr_solve::{Config, Solver};

/// Play a random game from 0^n for `choices.len()` rounds (or until no cell is
/// empty), using `choices` to pick the Pointer's cell and the Labeler's option.
fn random_reachable<V: Variant>(n: u32, choices: &[(u32, u32)]) -> Board<u64> {
    let s = Solver::<u64, V>::new(Config::default());
    let mut x = Board::<u64>::empty(n);
    for &(cj, co) in choices {
        let empt: Vec<u32> = spr_core::cells_of(x.empties()).collect();
        if empt.is_empty() {
            break;
        }
        let j = empt[cj as usize % empt.len()];
        let opts = s.labeler_options(x, j);
        x = opts[co as usize % opts.len()].2;
    }
    x
}

fn check_reachable_vs_brute<V: Variant>(n: u32, choices: &[(u32, u32)], r: u32) {
    let x = random_reachable::<V>(n, choices);
    let mut brute = Brute::new(V::KIND);
    let mut solver = Solver::<u64, V>::new(Config::default());
    assert_eq!(solver.value(x, r).unwrap(), brute.value(&x.to_cells(), r), "{} x={x} r={r}", V::KIND);
}

fn check_value_invariants<V: Variant>(n: u32, choices: &[(u32, u32)], r: u32) {
    let x = random_reachable::<V>(n, choices);
    let mut solver = Solver::<u64, V>::new(Config::default());
    let c = x.count();
    let v = solver.value(x, r).unwrap();
    let v_next = solver.value(x, r + 1).unwrap();
    assert!(c <= v && v <= c + r, "{} x={x} r={r} v={v}", V::KIND);
    assert!(v <= n);
    assert!(v <= v_next, "{} not monotone in r: x={x} r={r}", V::KIND);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(400))]

    #[test]
    fn reachable_matches_brute(n in 1u32..=8, choices in prop::collection::vec((0u32..64, 0u32..1024), 0..12), r in 0u32..=8) {
        check_reachable_vs_brute::<Full>(n, &choices, r.min(6));
        check_reachable_vs_brute::<Forced>(n, &choices, r);
        check_reachable_vs_brute::<Qv>(n, &choices, r);
    }

    #[test]
    fn value_invariants(n in 1u32..=14, choices in prop::collection::vec((0u32..64, 0u32..1024), 0..14), r in 0u32..=10) {
        check_value_invariants::<Forced>(n, &choices, r);
        if n <= 12 {
            check_value_invariants::<Qv>(n, &choices, r);
        }
        if n <= 11 {
            check_value_invariants::<Full>(n, &choices, r.min(8));
        }
    }

    #[test]
    fn opt_invariants(n in 1u32..=12, s in 1u32..=12) {
        // QV without its own reductions is the slowest variant at large s.
        let (n, s) = (n.min(12), s.min(12));
        let mut forced = Solver::<u64, Forced>::new(Config::default());
        let mut qv = Solver::<u64, Qv>::new(Config::default());
        let f = forced.opt(n, s).unwrap();
        let q = qv.opt(n, s).unwrap();
        prop_assert!(f <= n.min(s));
        prop_assert!(q <= f, "opt_qv > opt_forced at n={} s={}", n, s);
        prop_assert!(forced.opt(n, s + 1).unwrap() >= f);
        prop_assert!(forced.opt(n + 1, s).unwrap() >= f);
        prop_assert!(qv.opt(n, s + 1).unwrap() >= q);
        prop_assert!(qv.opt(n + 1, s).unwrap() >= q);
        if n <= 10 && s <= 9 {
            let mut full = Solver::<u64, Full>::new(Config::default());
            let v = full.opt(n, s).unwrap();
            prop_assert!(v <= n.min(s));
            prop_assert!(q <= v);
            prop_assert!(full.opt(n + 1, s).unwrap() >= v);
            prop_assert!(full.opt(n, s + 1).unwrap() >= v);
        }
    }
}

/// Monotonicity bounds on or off must give identical tables.
#[test]
fn monotone_flag_does_not_change_values() {
    fn run<V: Variant>(n_max: u32, s_max: u32) {
        let mut on = Solver::<u64, V>::new(Config { monotone: true, ..Config::default() });
        let mut off = Solver::<u64, V>::new(Config { monotone: false, ..Config::default() });
        for n in 1..=n_max {
            for s in 1..=s_max {
                assert_eq!(on.opt(n, s), off.opt(n, s), "{} n={n} s={s}", V::KIND);
            }
        }
    }
    run::<Forced>(16, 20);
    run::<Qv>(12, 12);
    run::<Full>(11, 12);
}

/// u64 and u128 boards give the same values.
#[test]
fn word_widths_agree() {
    let mut a = Solver::<u64, Forced>::new(Config::default());
    let mut b = Solver::<u128, Forced>::new(Config::default());
    for n in 1..=14 {
        for s in 1..=14 {
            assert_eq!(a.opt(n, s), b.opt(n, s));
        }
    }
}

/// (d) Qiao-Valiant Lemma 5: opt_qv(2^t - 1, t) = t for t <= 5.
#[test]
fn qiao_valiant_lemma_5() {
    for t in 1..=5u32 {
        let n = (1 << t) - 1;
        let mut qv = Solver::<u64, Qv>::new(Config::default());
        assert_eq!(qv.opt(n, t), Some(t), "opt_qv({n}, {t})");
    }
}
