//! The gap lemma V(x, r) >= |x| + opt(g(x), r), tested in the forced variant.
//!
//! Under forced erasure every reachable board has all '+' strictly left of all '-'.
//! The gap is the run of empty cells strictly between the last '+' and the first '-'
//! (from the left edge if there is no '+', to the right edge if there is no '-'); on such
//! a board every cell in that run is empty. The lemma claims
//! `V(x, r) >= |x| + opt(g(x), r)`, with `opt = opt_forced`. Nothing here assumes it:
//! every value is computed exactly by the forced solver, including `opt(g, r)` for
//! `r > 2g`, which the stored tables do not cover.

use crate::solver::Solver;
use crate::strategy::principal_variation;
use spr_core::{Board, Forced, MINUS, PLUS, Word, cells_of};
use std::collections::HashMap;

/// `(first cell of the gap, gap length)`, or `None` if some '+' is not left of every '-'.
pub fn gap(x: &Board<u64>) -> Option<(u32, u32)> {
    let plus = x.plus();
    let minus = x.minus();
    let last_plus: i64 = if plus == 0 { -1 } else { ((63 - plus.leading_zeros()) / 2) as i64 };
    let first_minus: i64 = if minus == 0 { x.len as i64 } else { (minus.trailing_zeros() / 2) as i64 };
    if last_plus >= first_minus {
        return None;
    }
    Some(((last_plus + 1) as u32, (first_minus - last_plus - 1) as u32))
}

/// A move from a board: `(parent bits, pointer, sign)`.
type Parent = (u64, u32, u8);

/// Every board reachable from `0^n` in the forced game, with one parent pointer each.
pub fn reachable(n: u32) -> HashMap<u64, Option<Parent>> {
    let mut seen: HashMap<u64, Option<Parent>> = HashMap::new();
    let start = Board::<u64>::empty(n);
    seen.insert(start.bits, None);
    let mut frontier = vec![start];
    while let Some(x) = frontier.pop() {
        for j in cells_of(x.empties()) {
            let e = x.erasable(j);
            for sign in [PLUS, MINUS] {
                let y = x.apply(j, e, sign, false);
                if let std::collections::hash_map::Entry::Vacant(v) = seen.entry(y.bits) {
                    v.insert(Some((x.bits, j, sign)));
                    frontier.push(y);
                }
            }
        }
    }
    seen
}

/// The move sequence from `0^n` to `bits`.
pub fn path_to(seen: &HashMap<u64, Option<Parent>>, n: u32, bits: u64) -> Vec<(Board<u64>, u32, u8)> {
    let mut out = Vec::new();
    let mut cur = bits;
    while let Some(Some((p, j, sign))) = seen.get(&cur) {
        out.push((Board { bits: *p, len: n }, *j, *sign));
        cur = *p;
    }
    out.reverse();
    out
}

pub struct Violation {
    pub n: u32,
    pub r: u32,
    pub board: Board<u64>,
    pub v: u32,
    pub g: u32,
    pub opt_g: u32,
    pub path: Vec<(Board<u64>, u32, u8)>,
}

pub struct CheckRow {
    pub n: u32,
    pub boards: usize,
    pub checks: u64,
    pub tight: u64,
    pub max_slack: u32,
}

/// Exhaustive check for one `n`: every reachable board `x` and every `r >= 1` with
/// `|x| + r <= 2n`, smallest `r` first. Returns the row statistics, or the first
/// violation. `None` if the solver hit its state limit.
pub fn check_n(solver: &mut Solver<u64, Forced>, n: u32) -> Option<Result<CheckRow, Violation>> {
    let seen = reachable(n);
    let mut boards: Vec<u64> = seen.keys().cloned().collect();
    boards.sort_unstable();
    let mut row = CheckRow { n, boards: boards.len(), checks: 0, tight: 0, max_slack: 0 };
    for r in 1..=2 * n {
        for &bits in &boards {
            let x = Board { bits, len: n };
            let c = x.count();
            if c + r > 2 * n {
                continue;
            }
            let (_, g) = gap(&x).expect("forced boards have '+' left of '-'");
            let v = solver.value(x, r)?;
            let og = solver.opt(g, r)?;
            row.checks += 1;
            if v < c + og {
                return Some(Err(Violation { n, r, board: x, v, g, opt_g: og, path: path_to(&seen, n, bits) }));
            }
            let slack = v - c - og;
            if slack == 0 {
                row.tight += 1;
            }
            row.max_slack = row.max_slack.max(slack);
        }
    }
    Some(Ok(row))
}

/// One position of a principal variation, with the gap decomposition.
pub struct PvRow {
    pub step: u32,
    pub r: u32,
    pub board: Board<u64>,
    pub v: u32,
    pub count: u32,
    pub g: u32,
    pub opt_g: u32,
    pub slack: u32,
    /// "stop", or where the Pointer points: "gap", "left" (left of the gap) or "right".
    pub region: &'static str,
    pub j: Option<u32>,
    pub refill: bool,
    pub erased: Vec<u32>,
    pub sign: Option<u8>,
    /// For the sign written at this step: `None` if it survives to the end, else the
    /// step at which it is erased.
    pub erased_at: Option<u32>,
}

/// A principal variation with fixed tie-breaks (centre-most optimal pointer, slowest-defeat
/// Labeler) from `0^n` with `s` rounds, annotated with the gap decomposition and the
/// fate of every sign. A refill is a pointer at a cell that an earlier erasure emptied
/// and that lies in the current gap.
pub fn pv_rows(solver: &mut Solver<u64, Forced>, n: u32, s: u32) -> Option<Vec<PvRow>> {
    let steps = principal_variation(solver, Board::<u64>::empty(n), s)?;
    let mut rows = Vec::new();
    let mut emptied_by_erasure = vec![false; n as usize];
    // cell -> index in `rows` of the step whose sign currently sits there
    let mut owner: Vec<Option<usize>> = vec![None; n as usize];
    for (k, st) in steps.iter().enumerate() {
        let x = st.board;
        let c = x.count();
        let (g0, g) = gap(&x).expect("forced boards have '+' left of '-'");
        let og = solver.opt(g, st.r)?;
        let slack = st.value - c - og;
        let mut row = PvRow {
            step: k as u32,
            r: st.r,
            board: x,
            v: st.value,
            count: c,
            g,
            opt_g: og,
            slack,
            region: "stop",
            j: None,
            refill: false,
            erased: Vec::new(),
            sign: None,
            erased_at: None,
        };
        if let Some((j, erased, sign)) = &st.action {
            let j = *j;
            row.region = if j < g0 {
                "left"
            } else if j < g0 + g {
                "gap"
            } else {
                "right"
            };
            row.refill = row.region == "gap" && emptied_by_erasure[j as usize];
            row.j = Some(j);
            row.erased = erased.clone();
            row.sign = Some(*sign);
            for &i in erased {
                if let Some(o) = owner[i as usize].take() {
                    let o: &mut PvRow = &mut rows[o];
                    o.erased_at = Some(k as u32);
                }
                emptied_by_erasure[i as usize] = true;
            }
            emptied_by_erasure[j as usize] = false;
            owner[j as usize] = Some(rows.len());
        }
        rows.push(row);
    }
    Some(rows)
}

/// Text form of a board, `0 + -`.
pub fn board_text<W: Word>(x: &Board<W>) -> String {
    x.to_string()
}
