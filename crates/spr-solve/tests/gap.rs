//! Sanity checks for the gap machinery.

use spr_core::{Board, Forced, VariantKind};
use spr_solve::brute::Brute;
use spr_solve::gap::{gap, reachable};
use spr_solve::{Config, Solver};

#[test]
fn gap_on_hand_boards() {
    let g = |s: &str| gap(&Board::<u64>::parse(s).unwrap());
    assert_eq!(g("00000"), Some((0, 5)));
    assert_eq!(g("+0-"), Some((1, 1)));
    assert_eq!(g("0+00-0"), Some((2, 2)));
    assert_eq!(g("-00"), Some((0, 0)));
    assert_eq!(g("00+"), Some((3, 0)));
    assert_eq!(g("+-"), Some((1, 0)));
    assert_eq!(g("-+"), None);
}

/// Every forced-reachable board has '+' left of '-', and the BFS finds all such boards
/// for small n (each one is reachable by writing its signs in order).
#[test]
fn reachable_boards_are_forced_shaped() {
    for n in 1..=8u32 {
        let seen = reachable(n);
        for &bits in seen.keys() {
            assert!(gap(&Board { bits, len: n }).is_some());
        }
        // boards with all '+' left of all '-': sum over split points of 2^n-type counts
        let shaped = (0..3u64.pow(n))
            .filter(|&k| {
                let mut v = Vec::new();
                let mut k = k;
                for _ in 0..n {
                    v.push((k % 3) as u8);
                    k /= 3;
                }
                Board::<u64>::from_cells(&v).plus_left_of_minus()
            })
            .count();
        assert_eq!(seen.len(), shaped, "n={n}");
    }
}

/// The lemma recomputed with the unreduced brute force for n <= 6 (independent of the
/// optimised solver), on every reachable board and every r with |x| + r <= 2n.
#[test]
fn lemma_on_brute_force_small_n() {
    let mut brute = Brute::new(VariantKind::Forced);
    let mut solver = Solver::<u64, Forced>::new(Config::default());
    for n in 1..=6u32 {
        for &bits in reachable(n).keys() {
            let x = Board::<u64> { bits, len: n };
            let (_, g) = gap(&x).unwrap();
            for r in 1..=(2 * n).saturating_sub(x.count()) {
                let v = brute.value(&x.to_cells(), r);
                let og = brute.opt(g as usize, r);
                assert_eq!(v, solver.value(x, r).unwrap());
                assert!(v >= x.count() + og, "violation n={n} x={x} r={r}");
            }
        }
    }
}
