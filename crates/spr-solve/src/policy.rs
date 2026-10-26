//! Exact value of a fixed Pointer policy against a best-responding Labeler.
//!
//! The Pointer's moves are fixed by a policy (a function of the board and the rounds
//! remaining); only its decision to stop is optimised. So
//!
//! ```text
//! W(x, 0) = |x|
//! W(x, r) = max( |x|, min over Labeler replies to policy(x, r) of W(x', r - 1) )
//! ```
//!
//! (just `|x|` if the policy has no move). Since stopping is always allowed and the
//! policy is one particular Pointer strategy, `W(0^n, s) <= opt(n, s)` in the same
//! variant. Against a Labeler that may erase any subset of the erasable set (`full`),
//! `W` is a lower bound on `opt_full`, with no dependence on the erasure observation.
//!
//! Memoised on the raw board (no reflection or banking, since a policy need not be
//! invariant under either).

use hashbrown::HashMap;
use spr_core::{Board, MINUS, PLUS, cells_of};

/// A Pointer policy: the cell to point at, or `None` to stop. It must return an empty
/// cell. It is evaluated on every board the Labeler can produce, including boards
/// where some '+' lies right of some '-' (possible against a full-game Labeler).
pub trait Policy {
    fn name(&self) -> String;
    fn choose(&self, x: &Board<u64>, r: u32) -> Option<u32>;
}

/// Last '+' and first '-' positions (`-1` and `n` when absent).
pub fn last_plus_first_minus(x: &Board<u64>) -> (i64, i64) {
    let plus = x.plus();
    let minus = x.minus();
    let lp = if plus == 0 { -1 } else { ((63 - plus.leading_zeros()) / 2) as i64 };
    let fm = if minus == 0 { x.len as i64 } else { (minus.trailing_zeros() / 2) as i64 };
    (lp, fm)
}

/// Empty cells strictly between `a` and `b`.
fn empties_between(x: &Board<u64>, a: i64, b: i64) -> Vec<u32> {
    cells_of(x.empties()).filter(|&c| (c as i64) > a && (c as i64) < b).collect()
}

/// Binary search: point at the middle of the safe interval (the empty cells between the
/// last '+' and the first '-', which exists only when every '+' is left of every '-');
/// stop when it is empty. Every move is to a cell whose erasable set is empty, so a
/// full-game Labeler has no erasure choice and this equals the forced game.
pub struct BinarySearch;

impl Policy for BinarySearch {
    fn name(&self) -> String {
        "binary-search".into()
    }
    fn choose(&self, x: &Board<u64>, _r: u32) -> Option<u32> {
        let (lp, fm) = last_plus_first_minus(x);
        if lp >= fm {
            return None;
        }
        let g = empties_between(x, lp, fm);
        if g.is_empty() {
            return None;
        }
        // All cells between lp and fm are empty; take the middle one.
        Some(g[(g.len() - 1) / 2])
    }
}

/// Where to point inside a non-empty gap.
#[derive(Clone, Copy, Debug)]
pub enum GapRule {
    /// Next to the first '-' (refill from the right).
    Right,
    /// The middle of the gap.
    Centre,
    /// Next to the last '+'.
    Left,
}

/// The walker read off the long-game principal variations: refill the gap; when there
/// is no gap, jump to the first empty cell right of the first '-', which (under forced
/// erasure) erases the '-' block and reopens it as a gap. Defined on every board: when
/// some '+' lies right of some '-', there is no gap and the jump rule applies.
pub struct Walker {
    pub gap: GapRule,
}

impl Policy for Walker {
    fn name(&self) -> String {
        format!("walker-{:?}", self.gap).to_lowercase()
    }
    fn choose(&self, x: &Board<u64>, _r: u32) -> Option<u32> {
        let (lp, fm) = last_plus_first_minus(x);
        if lp < fm {
            let g = empties_between(x, lp, fm);
            if !g.is_empty() {
                return Some(match self.gap {
                    GapRule::Right => *g.last().unwrap(),
                    GapRule::Left => g[0],
                    GapRule::Centre => g[(g.len() - 1) / 2],
                });
            }
        }
        // Jump: first empty cell right of the first '-'.
        cells_of(x.empties()).find(|&c| (c as i64) > fm)
    }
}

/// The walker with fresh starts: in a gap with no '-' to its right (a fresh region),
/// point at its centre; in a gap bounded on the right by a '-' block, refill next to that
/// block; with no gap, jump to the first empty cell right of the first '-'. This is the
/// rightward walk seen in the long-game principal variations (the walking block of the long lines).
pub struct WalkerFresh;

impl Policy for WalkerFresh {
    fn name(&self) -> String {
        "walker-fresh".into()
    }
    fn choose(&self, x: &Board<u64>, _r: u32) -> Option<u32> {
        let n = x.len as i64;
        let (lp, fm) = last_plus_first_minus(x);
        if lp < fm {
            let g = empties_between(x, lp, fm);
            if !g.is_empty() {
                return Some(if fm == n { g[(g.len() - 1) / 2] } else { *g.last().unwrap() });
            }
        }
        cells_of(x.empties()).find(|&c| (c as i64) > fm)
    }
}

/// Two-sided walker: walks towards the side with more free room. With `lp` the last '+'
/// and `fm` the first '-', the room on the right is `n - 1 - fm` and on the left is `lp`.
/// On the empty board it points at the centre. Walking right, it refills the gap next to
/// the '-' block and, with no gap, jumps to the first empty cell right of `fm`. Walking
/// left is the mirror image: refill next to the '+' block, jump to the last empty cell
/// left of `lp`. If the preferred jump has no empty cell it tries the other side, then
/// stops.
pub struct WalkerBoth;

impl Policy for WalkerBoth {
    fn name(&self) -> String {
        "walker-both".into()
    }
    fn choose(&self, x: &Board<u64>, _r: u32) -> Option<u32> {
        let n = x.len as i64;
        let (lp, fm) = last_plus_first_minus(x);
        let right = (n - 1 - fm) >= lp;
        if lp < fm {
            let g = empties_between(x, lp, fm);
            if !g.is_empty() {
                if lp == -1 && fm == n {
                    return Some(g[(g.len() - 1) / 2]);
                }
                return Some(if right { *g.last().unwrap() } else { g[0] });
            }
        }
        let jr = cells_of(x.empties()).find(|&c| (c as i64) > fm);
        let jl = cells_of(x.empties()).filter(|&c| (c as i64) < lp).last();
        if right { jr.or(jl) } else { jl.or(jr) }
    }
}

pub struct Evaluator<'a> {
    policy: &'a dyn Policy,
    full: bool,
    memo: HashMap<(u64, u32, u32), u8>,
    max_states: usize,
    pub nodes: u64,
    aborted: bool,
}

impl<'a> Evaluator<'a> {
    /// `full`: the Labeler may erase any subset; otherwise it must erase all.
    pub fn new(policy: &'a dyn Policy, full: bool, max_states: usize) -> Self {
        Evaluator { policy, full, memo: HashMap::new(), max_states, nodes: 0, aborted: false }
    }

    pub fn states(&self) -> usize {
        self.memo.len()
    }

    /// `W(0^n, s)`, or `None` if the state limit was reached.
    pub fn value_root(&mut self, n: u32, s: u32) -> Option<u32> {
        let v = self.value(Board::empty(n), s);
        (!self.aborted).then_some(v)
    }

    pub fn value(&mut self, x: Board<u64>, r: u32) -> u32 {
        if self.aborted {
            return 0;
        }
        let c = x.count();
        if r == 0 {
            return c;
        }
        let key = (x.bits, x.len, r);
        if let Some(&v) = self.memo.get(&key) {
            return v as u32;
        }
        self.nodes += 1;
        let mut best = c;
        if let Some(j) = self.policy.choose(&x, r) {
            debug_assert_eq!(x.cell(j), 0);
            let m = x.erasable(j);
            let mut worst = u32::MAX;
            let mut sub = m;
            'outer: loop {
                for sign in [PLUS, MINUS] {
                    let v = self.value(x.apply(j, sub, sign, false), r - 1);
                    if v < worst {
                        worst = v;
                        if worst <= best {
                            break 'outer;
                        }
                    }
                }
                if !self.full || sub == 0 {
                    break;
                }
                sub = sub.wrapping_sub(1) & m;
            }
            best = best.max(worst);
        }
        if self.memo.len() >= self.max_states {
            self.aborted = true;
            return best;
        }
        self.memo.insert(key, best as u8);
        best
    }
}

/// A Labeler reply chosen by `Evaluator` (for line extraction): the first reply that
/// attains the minimum.
pub fn worst_reply(ev: &mut Evaluator, x: Board<u64>, r: u32, j: u32) -> (u64, u8, Board<u64>, u32) {
    let m = x.erasable(j);
    let mut best: Option<(u64, u8, Board<u64>, u32)> = None;
    let mut sub = m;
    loop {
        for sign in [PLUS, MINUS] {
            let y = x.apply(j, sub, sign, false);
            let v = ev.value(y, r - 1);
            if best.as_ref().is_none_or(|b| v < b.3) {
                best = Some((sub, sign, y, v));
            }
        }
        if !ev.full || sub == 0 {
            break;
        }
        sub = sub.wrapping_sub(1) & m;
    }
    best.unwrap()
}
