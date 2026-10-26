//! Validation (b) and the reduction proofs, checked exhaustively against the
//! unreduced brute force.

use spr_core::naive;
use spr_core::{Board, Forced, Full, Qv, Variant, VariantKind, Word};
use spr_solve::brute::Brute;
use spr_solve::{Config, Solver};

/// Every board of length n over the given number of codes (3: no burned, 4: with).
fn all_boards(n: usize, codes: u8) -> impl Iterator<Item = Vec<u8>> {
    let total = (codes as u64).pow(n as u32);
    (0..total).map(move |mut k| {
        let mut v = vec![0u8; n];
        for c in v.iter_mut() {
            *c = (k % codes as u64) as u8;
            k /= codes as u64;
        }
        v
    })
}

fn codes(kind: VariantKind) -> u8 {
    if kind.burn() { 4 } else { 3 }
}

/// Largest r checked on all states. r = 10 exceeds n for every n <= 9, so each
/// board is checked well past the point where its value can still change.
fn r_max(_kind: VariantKind, _n: usize) -> u32 {
    10
}

/// Optimised solver equals brute force on every state (x, r), n <= 9.
fn solver_matches_brute_all_states<W: Word, V: Variant>(monotone: bool) {
    solver_matches_brute_with::<W, V>(Config { monotone, ..Config::default() });
}

fn solver_matches_brute_with<W: Word, V: Variant>(cfg: Config) {
    let monotone = cfg.monotone;
    let kind = V::KIND;
    for n in 0..=9usize {
        let mut brute = Brute::new(kind);
        let mut solver = Solver::<W, V>::new(cfg);
        let mut checked = 0u64;
        // Increasing r so the monotone bounds are actually exercised.
        for r in 0..=r_max(kind, n) {
            for x in all_boards(n, codes(kind)) {
                let b = brute.value(&x, r);
                let s = solver.value(Board::<W>::from_cells(&x), r).unwrap();
                assert_eq!(s, b, "{kind} x={x:?} r={r} monotone={monotone}");
                checked += 1;
            }
        }
        assert!(checked > 0);
    }
}

#[test]
fn full_all_states() {
    solver_matches_brute_all_states::<u64, Full>(true);
    solver_matches_brute_all_states::<u64, Full>(false);
}

#[test]
fn forced_all_states() {
    solver_matches_brute_all_states::<u64, Forced>(true);
    solver_matches_brute_all_states::<u64, Forced>(false);
    solver_matches_brute_all_states::<u128, Forced>(true);
}

#[test]
fn qv_all_states() {
    solver_matches_brute_all_states::<u64, Qv>(true);
    solver_matches_brute_all_states::<u64, Qv>(false);
    solver_matches_brute_with::<u64, Qv>(Config { qv_reduce: false, ..Config::default() });
    solver_matches_brute_with::<u128, Qv>(Config::default());
}

/// QV reduction identity on the brute force: V(x, r) = bank + V(reduced x, r).
#[test]
fn qv_reduction_holds_on_brute_force() {
    let mut brute = Brute::new(VariantKind::Qv);
    for n in 0..=9usize {
        for r in 0..=10 {
            for x in all_boards(n, 4) {
                let (bank, core) = naive::qv_reduce(&x);
                assert_eq!(bank + brute.value(&core, r), brute.value(&x, r), "qv reduce {x:?} r={r}");
            }
        }
    }
}

/// Reflection and banking identities checked on the brute force itself (no packed
/// code involved): V(x, r) = V(reflect x, r) and V(x, r) = p + q + V(core, r).
#[test]
fn reductions_hold_on_brute_force() {
    for kind in VariantKind::ALL {
        let mut brute = Brute::new(kind);
        for n in 0..=9usize {
            for r in 0..=r_max(kind, n) {
                for x in all_boards(n, codes(kind)) {
                    let v = brute.value(&x, r);
                    assert_eq!(brute.value(&naive::reflect(&x), r), v, "{kind} reflect {x:?} r={r}");
                    let (p, q) = naive::bank_lengths(&x);
                    let core = x[p..n - q].to_vec();
                    assert_eq!((p + q) as u32 + brute.value(&core, r), v, "{kind} bank {x:?} r={r}");
                }
            }
        }
    }
}

/// opt(n, s) from the empty board, n <= 8, s <= 2n, every variant.
#[test]
fn opt_matches_brute_n_le_8() {
    fn run<V: Variant>(s_cap: u32) {
        for n in 1..=8u32 {
            let mut brute = Brute::new(V::KIND);
            let mut solver = Solver::<u64, V>::new(Config::default());
            for s in 1..=(2 * n).min(s_cap) {
                assert_eq!(solver.opt(n, s), Some(brute.opt(n as usize, s)), "{} n={n} s={s}", V::KIND);
            }
        }
    }
    run::<Forced>(16);
    run::<Qv>(16);
    run::<Full>(16);
}
