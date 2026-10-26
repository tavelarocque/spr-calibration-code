//! Strategy extraction: every optimal Pointer move at a state, one principal
//! variation, and the reuse analysis of a forced-variant line. All values come from
//! the exact solver; no cutoffs are used to decide optimality.

use crate::solver::Solver;
use spr_core::{Board, Forced, Qv, Variant, Word, cells_of};

#[derive(Clone, Debug)]
pub struct RootAnalysis {
    pub value: u32,
    /// Stopping now scores `|x| = value`.
    pub stop_optimal: bool,
    /// `(j, Labeler value of j)` for every empty cell `j`.
    pub pointer_values: Vec<(u32, u32)>,
}

impl RootAnalysis {
    pub fn optimal_pointers(&self) -> Vec<u32> {
        self.pointer_values.iter().filter(|&&(_, v)| v == self.value).map(|&(j, _)| j).collect()
    }
}

pub fn analyse<W: Word, V: Variant>(solver: &mut Solver<W, V>, x: Board<W>, r: u32) -> Option<RootAnalysis> {
    let value = solver.value(x, r)?;
    let mut pointer_values = Vec::new();
    if r >= 1 {
        for j in cells_of(x.empties()) {
            pointer_values.push((j, solver.pointer_value(x, r, j)?));
        }
    }
    Some(RootAnalysis { value, stop_optimal: x.count() == value, pointer_values })
}

#[derive(Clone, Debug)]
pub struct Step<W: Word> {
    pub board: Board<W>,
    pub r: u32,
    pub value: u32,
    /// Every optimal pointer at this state (empty if `r = 0`).
    pub optimal: Vec<u32>,
    /// `None`: the Pointer stops. Otherwise `(j, erased cells, sign)`.
    pub action: Option<(u32, Vec<u32>, u8)>,
}

/// Least `t <= r` with `V(x, t) >= target`: how many rounds the Pointer needs from
/// `x` to secure `target` (`V` is non-decreasing in `t`). `None` if even `r` rounds
/// do not suffice, or if the table limit was hit.
pub fn rounds_needed<W: Word, V: Variant>(solver: &mut Solver<W, V>, x: Board<W>, r: u32, target: u32) -> Option<u32> {
    (0..=r).find(|&t| solver.value(x, t).is_some_and(|v| v >= target))
}

/// One principal variation. The Pointer stops as soon as stopping is optimal and
/// otherwise takes the optimal cell closest to the centre of the board (the lower
/// one on a tie). The Labeler plays the slowest defeat: among its optimal replies
/// it takes one that maximises the number of rounds the Pointer still needs to
/// secure the value, breaking ties in the order maximal erasure first, '+' before
/// '-'. Without this, ties let the Labeler hand out signs and the line stops early.
pub fn principal_variation<W: Word, V: Variant>(
    solver: &mut Solver<W, V>,
    x0: Board<W>,
    r0: u32,
) -> Option<Vec<Step<W>>> {
    let mut steps = Vec::new();
    let (mut x, mut r) = (x0, r0);
    loop {
        let a = analyse(solver, x, r)?;
        let optimal = a.optimal_pointers();
        if a.stop_optimal || r == 0 {
            steps.push(Step { board: x, r, value: a.value, optimal, action: None });
            return Some(steps);
        }
        let centre2 = x.len as i32 - 1;
        let j = *optimal.iter().min_by_key(|&&j| ((2 * j as i32 - centre2).abs(), j)).unwrap();
        let mut chosen: Option<(u32, W, u8, Board<W>)> = None;
        for (e, sign, child) in solver.labeler_options(x, j) {
            if solver.value(child, r - 1)? != a.value {
                continue;
            }
            let need = rounds_needed(solver, child, r - 1, a.value)?;
            if chosen.as_ref().is_none_or(|c| need > c.0) {
                chosen = Some((need, e, sign, child));
            }
        }
        let (_, e, sign, child) = chosen.expect("some Labeler option attains the value");
        steps.push(Step { board: x, r, value: a.value, optimal, action: Some((j, cells_of(e).collect(), sign)) });
        x = child;
        r -= 1;
    }
}

/// One step of a reuse line.
#[derive(Clone, Debug)]
pub struct ReuseStep {
    pub board: Board<u64>,
    pub r: u32,
    /// Spread mask of cells pointed at so far.
    pub used: u64,
    pub v_forced: u32,
    /// QV value of the same history: used-and-now-empty cells are burned.
    pub v_qv: u32,
    /// `(j, forced Labeler value of j)` for every empty `j`.
    pub pointer_values: Vec<(u32, u32)>,
    pub action: Option<(u32, Vec<u32>, u8)>,
}

impl ReuseStep {
    /// Best forced Labeler value over empty cells never pointed at, and over
    /// empty cells pointed at before (reused). `None` if there is no such cell.
    pub fn best_fresh_and_reused(&self) -> (Option<u32>, Option<u32>) {
        let pick = |reused: bool| {
            self.pointer_values.iter().filter(|&&(j, _)| (self.used >> (2 * j) & 1 == 1) == reused).map(|&(_, v)| v).max()
        };
        (pick(false), pick(true))
    }
}

/// The QV image of a forced-variant history: every cell that has been pointed at
/// and is now empty becomes burned. QV from this board is exactly the forced game
/// continued under the rule that a cell can be pointed at only once.
pub fn qv_image(x: Board<u64>, used: u64) -> Board<u64> {
    let m = used & x.empties();
    Board { bits: x.bits | m | (m << 1), len: x.len }
}

/// A forced-variant optimal line built to expose reuse. The Pointer plays an
/// optimal move, preferring cells never pointed at (lowest index first). The
/// Labeler plays an optimal reply, breaking ties towards the reply that makes the
/// QV image of the history worst for the Pointer. Along this line the Pointer
/// points at a reused cell only when no fresh cell is optimal.
pub fn reuse_line(
    forced: &mut Solver<u64, Forced>,
    qv: &mut Solver<u64, Qv>,
    n: u32,
    s: u32,
) -> Option<Vec<ReuseStep>> {
    let mut steps = Vec::new();
    let (mut x, mut r, mut used) = (Board::<u64>::empty(n), s, 0u64);
    loop {
        let a = analyse(forced, x, r)?;
        let v_qv = qv.value(qv_image(x, used), r)?;
        let mut step = ReuseStep {
            board: x,
            r,
            used,
            v_forced: a.value,
            v_qv,
            pointer_values: a.pointer_values.clone(),
            action: None,
        };
        if a.stop_optimal || r == 0 {
            steps.push(step);
            return Some(steps);
        }
        let opt = a.optimal_pointers();
        let j = *opt.iter().find(|&&j| used >> (2 * j) & 1 == 0).unwrap_or(&opt[0]);
        let used_next = used | (1u64 << (2 * j));
        let mut best: Option<(u32, u64, u8, Board<u64>)> = None;
        for (e, sign, child) in forced.labeler_options(x, j) {
            if forced.value(child, r - 1)? != a.value {
                continue;
            }
            let q = qv.value(qv_image(child, used_next), r - 1)?;
            if best.as_ref().is_none_or(|b| q < b.0) {
                best = Some((q, e, sign, child));
            }
        }
        let (_, e, sign, child) = best.expect("some Labeler option attains the value");
        step.action = Some((j, cells_of(e).collect(), sign));
        steps.push(step);
        x = child;
        used = used_next;
        r -= 1;
    }
}
