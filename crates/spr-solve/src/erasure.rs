//! The erasure lemma: is maximal erasure always optimal for the Labeler?
//!
//! For every state `(x, r)` reachable in the full game from `0^n` within `s_max`
//! rounds (any Pointer move, any Labeler response), and every empty `j`, compute
//! two independent exact quantities from exact child values, with no cutoffs:
//!
//! * `best_any = min over all (E, sign) of V_full(x', r - 1)`
//! * `best_max = min over sign of V_full(x' with E maximal, r - 1)`
//!
//! A counterexample is a state and pointer with `best_any < best_max`: some
//! non-maximal erasure is strictly better for the Labeler.
//!
//! If there is none on the whole reachable set, then `V_full = V_forced` on every
//! reachable state, by induction on `r`: with no strict improvement,
//! `V_full(x, r) = max(|x|, max_j min_sign V_full(x_max, r - 1))`, the maximal-erasure
//! children are themselves reachable at `r - 1`, and the forced recursion is the
//! same formula.
//!
//! The check at `(x, r)` depends only on the canonical banked core of `x` (both
//! quantities shift by the bank and are reflection invariant), so each distinct
//! `(core, r)` is checked once, across all `n`.

use crate::solver::Solver;
use hashbrown::{HashMap, HashSet};
use spr_core::{Board, Full, Word, cells_of};

#[derive(Clone, Debug)]
pub struct Move<W: Word> {
    pub before: Board<W>,
    pub j: u32,
    pub erased: Vec<u32>,
    pub sign: u8,
}

#[derive(Clone, Debug)]
pub struct Counterexample<W: Word> {
    pub n: u32,
    pub r: u32,
    pub board: Board<W>,
    pub j: u32,
    pub best_any: u32,
    pub best_any_option: (Vec<u32>, u8),
    pub best_max: u32,
    /// Moves from `0^n` reaching `board` (shortest).
    pub path: Vec<Move<W>>,
}

#[derive(Clone, Debug, Default)]
pub struct RowStats {
    pub n: u32,
    pub s_max: u32,
    /// Distinct reachable boards (raw, not reduced) at depth <= s_max - 1.
    pub boards: u64,
    /// Reachable `(x, r)` pairs with `r >= 1` covered.
    pub pairs: u64,
    /// Distinct `(core, r)` pairs actually checked in this row (new ones only).
    pub core_checks: u64,
    /// `(state, pointer)` comparisons performed in this row.
    pub pointer_checks: u64,
}

pub struct ErasureSearch<W: Word> {
    pub solver: Solver<W, Full>,
    checked: HashSet<(W, u32, u32)>,
}

type Parent<W> = (u32, W, u32, W, u8);

impl<W: Word> ErasureSearch<W> {
    pub fn new(solver: Solver<W, Full>) -> Self {
        ErasureSearch { solver, checked: HashSet::new() }
    }

    /// Check every reachable state of `SPR(n, s)` for all `s <= s_max`. Returns the
    /// row statistics and the counterexample with smallest `r` in this row, if any.
    /// `None` if the solver ran out of table space.
    pub fn check_row(&mut self, n: u32, s_max: u32) -> Option<(RowStats, Option<Counterexample<W>>)> {
        // Breadth-first search over raw boards with parent pointers:
        // bits -> (depth, parent bits, j, erased mask, sign).
        let start = Board::<W>::empty(n);
        let mut parent: HashMap<W, Parent<W>> = HashMap::new();
        parent.insert(start.bits, (0, start.bits, 0, W::ZERO, 0));
        let mut frontier = vec![start];
        let mut by_depth: Vec<Vec<Board<W>>> = vec![vec![start]];
        for d in 1..s_max {
            let mut next = Vec::new();
            for &x in &frontier {
                for j in cells_of(x.empties()) {
                    for (e, sign, y) in self.solver.labeler_options(x, j) {
                        if !parent.contains_key(&y.bits) {
                            parent.insert(y.bits, (d, x.bits, j, e, sign));
                            next.push(y);
                        }
                    }
                }
            }
            if next.is_empty() {
                break;
            }
            by_depth.push(next.clone());
            frontier = next;
        }

        let mut stats = RowStats { n, s_max, boards: parent.len() as u64, ..RowStats::default() };
        // Smallest r first, so the first counterexample found has minimal r.
        for r in 1..=s_max {
            for (d, boards) in by_depth.iter().enumerate() {
                if d as u32 + r > s_max {
                    break;
                }
                for &x in boards {
                    stats.pairs += 1;
                    if x.empties().is_zero() {
                        continue;
                    }
                    let (_, core) = x.reduce();
                    if !self.checked.insert((core.bits, core.len, r)) {
                        continue;
                    }
                    stats.core_checks += 1;
                    if let Some(mut cx) = self.check_state(x, r, &mut stats)? {
                        cx.n = n;
                        cx.path = Self::path(&parent, start, x);
                        return Some((stats, Some(cx)));
                    }
                }
            }
        }
        Some((stats, None))
    }

    /// Compare best-subset and maximal-erasure values at `(x, r)` for every `j`.
    fn check_state(&mut self, x: Board<W>, r: u32, stats: &mut RowStats) -> Option<Option<Counterexample<W>>> {
        for j in cells_of(x.empties()) {
            stats.pointer_checks += 1;
            let m = x.erasable(j);
            let mut best_any = u32::MAX;
            let mut best_any_option = (Vec::new(), 0);
            let mut best_max = u32::MAX;
            for (e, sign, child) in self.solver.labeler_options(x, j) {
                let v = self.solver.value(child, r - 1)?;
                if v < best_any {
                    best_any = v;
                    best_any_option = (cells_of(e).collect(), sign);
                }
                if e == m {
                    best_max = best_max.min(v);
                }
            }
            if best_any < best_max {
                return Some(Some(Counterexample {
                    n: x.len,
                    r,
                    board: x,
                    j,
                    best_any,
                    best_any_option,
                    best_max,
                    path: Vec::new(),
                }));
            }
        }
        Some(None)
    }

    fn path(parent: &HashMap<W, Parent<W>>, start: Board<W>, x: Board<W>) -> Vec<Move<W>> {
        let mut moves = Vec::new();
        let mut cur = x.bits;
        while cur != start.bits {
            let &(_, p, j, e, sign) = &parent[&cur];
            moves.push(Move { before: Board { bits: p, len: x.len }, j, erased: cells_of(e).collect(), sign });
            cur = p;
        }
        moves.reverse();
        moves
    }
}
