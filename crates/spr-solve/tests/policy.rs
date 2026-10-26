//! The fixed-policy evaluator.

use spr_core::{Forced, Full};
use spr_solve::policy::{BinarySearch, Evaluator, GapRule, Policy, Walker, WalkerBoth, WalkerFresh};
use spr_solve::{Config, Solver};

fn r_bs(n: u32, s: u32) -> u32 {
    if n == 0 || s == 0 {
        return 0;
    }
    r_bs(n, s - 1).max(1 + r_bs((n - 1) / 2, s - 1))
}

/// Binary search points only at safe cells, so the full-game Labeler has no erasure
/// choice: its value is R(n, s) in both variants.
#[test]
fn binary_search_is_r_in_both_variants() {
    for full in [false, true] {
        for n in 1..=12 {
            for s in 1..=2 * n {
                let mut ev = Evaluator::new(&BinarySearch, full, 10_000_000);
                assert_eq!(ev.value_root(n, s), Some(r_bs(n, s)), "n={n} s={s} full={full}");
            }
        }
    }
}

/// A fixed policy is one Pointer strategy: its value never exceeds opt in its variant,
/// and the full-game Labeler can only do better than the forced one.
#[test]
fn policies_never_beat_opt() {
    let policies: Vec<Box<dyn Policy>> = vec![
        Box::new(Walker { gap: GapRule::Right }),
        Box::new(Walker { gap: GapRule::Centre }),
        Box::new(Walker { gap: GapRule::Left }),
        Box::new(WalkerFresh),
        Box::new(WalkerBoth),
    ];
    let mut forced = Solver::<u64, Forced>::new(Config::default());
    let mut full = Solver::<u64, Full>::new(Config::default());
    for p in &policies {
        for n in 1..=9 {
            for s in 1..=2 * n {
                let wf = Evaluator::new(p.as_ref(), false, 10_000_000).value_root(n, s).unwrap();
                let wu = Evaluator::new(p.as_ref(), true, 10_000_000).value_root(n, s).unwrap();
                assert!(wf <= forced.opt(n, s).unwrap(), "{} forced n={n} s={s}", p.name());
                assert!(wu <= full.opt(n, s).unwrap(), "{} full n={n} s={s}", p.name());
                assert!(wu <= wf, "{} n={n} s={s}", p.name());
            }
        }
    }
}
