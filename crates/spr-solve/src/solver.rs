//! Memoised depth-first minimax with a transposition table keyed by
//! `(canonical banked core, r)`.
//!
//! # Exactness
//!
//! Every value written to the table is the exact `V(core, r)`. The cutoffs are
//! local to one Pointer node and never narrow the window of a child search:
//!
//! * Pointer node: `best` starts at `|x|` (the stop value) and the loop over `j`
//!   ends once `best` reaches an upper bound `ub` on `V(x, r)`. Since
//!   `best <= V <= ub`, the result is then exact.
//! * Labeler node for pointer `j`: the loop over `(E, sign)` ends as soon as one
//!   child is `<= best`. That `j` then cannot raise `best`, which is all the parent
//!   needs. The partial minimum is never stored, because Labeler nodes are not
//!   memoised.
//!
//! The upper bound is `ub = |x| + min(r, #empty)`. Each round writes one sign into
//! a cell that was empty and erases nothing else into a sign, so `|x|` grows by at
//! most one per round, and it can only grow by pointing at an empty cell. When every
//! empty cell is a fresh one (full and forced), this is `min(|x| + r, n)` as in
//! the task. In QV, `#empty` counts fresh cells only, which is tighter but just as
//! valid: a burned cell never holds a sign again.
//!
//! The order in which the Pointer's cells are tried (centre first) does not affect
//! the result: the node returns `max(|x|, max_j L_j)` whatever the order, because
//! each `L_j` is either computed exactly or shown to be `<= best`, and `best` is
//! at every moment a lower bound on that maximum. The order only changes how
//! early `best` rises, and so how many Labeler cutoffs fire.
//!
//! Before searching any child, the Labeler makes a cheap pass over its options
//! (`labeler`). If some child has an upper bound `<= best` that needs no
//! search (the trivial bound, or a stored value), the minimum is `<= best` and the
//! pointer is refuted at once.
//!
//! # Monotonicity in r (optional, `Config::monotone`)
//!
//! `V(x, r) >= V(x, r - 1)`: with `r` rounds the Pointer can play an optimal
//! `(r - 1)`-round strategy and then stop. So a stored `V(x, r - 1)` is a lower
//! bound and a stored `V(x, r + 1)` an upper bound on `V(x, r)`. They are used only
//! to prune: the lower bound raises the starting `best` (which tightens the Labeler
//! cutoff), the upper bound lowers `ub`. The node's value is still computed as
//! `max(lb, |x|, max_j L_j)` with `L_j` the exact Labeler values above `best`, which
//! equals `V(x, r)` because `V >= lb` and `V = max(|x|, max_j L_j)`. A bound is
//! never written to the table as if it were a value at another `r`.
//!
//! # Running out of memory
//!
//! When the table reaches `Config::max_states` the solver sets a sticky `aborted`
//! flag. From then on it stores nothing, and every search returns a meaningless
//! number that the public entry points turn into `None`. Values stored before the
//! abort are exact and stay usable.

use hashbrown::HashTable;
use spr_core::{Board, MINUS, PLUS, Variant, Word, cells_of};
use std::marker::PhantomData;

#[derive(Clone, Copy, Debug)]
pub struct Config {
    /// Use stored values at r - 1 and r + 1 as bounds.
    pub monotone: bool,
    /// Stop storing (and abort) once the table holds this many entries.
    pub max_states: usize,
    /// Allocate the whole table up front, avoiding the 3x peak of a resize.
    pub reserve: bool,
    /// QV only: delete burned cells and bank permanent signs
    /// ([`Board::reduce_qv`]) instead of plain banking.
    pub qv_reduce: bool,
}

impl Default for Config {
    fn default() -> Self {
        Config { monotone: true, max_states: usize::MAX, reserve: false, qv_reduce: true }
    }
}

/// One table entry. `packed(4)` keeps the u64 entry at 12 bytes.
#[derive(Clone, Copy)]
#[repr(C, packed(4))]
struct Entry<W: Word> {
    bits: W,
    len: u8,
    r: u8,
    val: u8,
}

#[inline(always)]
fn fmix64(mut k: u64) -> u64 {
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51_afd7_ed55_8ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    k ^ (k >> 33)
}

#[inline(always)]
fn hash_key<W: Word>(bits: W, len: u32, r: u32) -> u64 {
    let (lo, hi) = bits.halves();
    fmix64(lo ^ fmix64(hi ^ ((len as u64) << 8 | r as u64)))
}

pub struct Solver<W: Word, V: Variant> {
    table: HashTable<Entry<W>>,
    cfg: Config,
    /// Pointer-node visits (calls to the internal search), cumulative.
    pub nodes: u64,
    aborted: bool,
    _v: PhantomData<V>,
}

impl<W: Word, V: Variant> Solver<W, V> {
    pub fn new(cfg: Config) -> Self {
        let table = if cfg.reserve && cfg.max_states != usize::MAX {
            HashTable::with_capacity(cfg.max_states)
        } else {
            HashTable::new()
        };
        Solver { table, cfg, nodes: 0, aborted: false, _v: PhantomData }
    }

    pub fn states(&self) -> usize {
        self.table.len()
    }

    pub fn table_bytes(&self) -> usize {
        self.table.allocation_size()
    }

    /// Number of stored entries per `r`, and per core length.
    pub fn histogram(&self) -> (Vec<u64>, Vec<u64>) {
        let (mut by_r, mut by_len) = (vec![0u64; 256], vec![0u64; 65]);
        for e in self.table.iter() {
            let (r, len) = (e.r, e.len);
            by_r[r as usize] += 1;
            by_len[len as usize] += 1;
        }
        while by_r.last() == Some(&0) {
            by_r.pop();
        }
        while by_len.last() == Some(&0) {
            by_len.pop();
        }
        (by_r, by_len)
    }

    pub fn aborted(&self) -> bool {
        self.aborted
    }

    /// Exact `V(x, r)` for any board, or `None` if the table limit was hit.
    pub fn value(&mut self, x: Board<W>, r: u32) -> Option<u32> {
        let (bank, core) = self.reduce(x);
        let v = self.search(core, r);
        (!self.aborted).then_some(bank + v)
    }

    /// `opt(n, s) = V(0^n, s)`.
    pub fn opt(&mut self, n: u32, s: u32) -> Option<u32> {
        self.value(Board::empty(n), s)
    }

    /// Exact Labeler value of pointer `j` at `(x, r)`, `r >= 1`: the minimum over
    /// every allowed `(E, sign)` of `V(x', r - 1)`. No cutoff.
    pub fn pointer_value(&mut self, x: Board<W>, r: u32, j: u32) -> Option<u32> {
        let mut worst = u32::MAX;
        for (_, _, child) in self.labeler_options(x, j) {
            worst = worst.min(self.value(child, r - 1)?);
        }
        Some(worst)
    }

    /// Every Labeler option `(E, sign, x')` for pointer `j`, maximal erasure first.
    /// The empty erasure set is included in the full variant.
    pub fn labeler_options(&self, x: Board<W>, j: u32) -> Vec<(W, u8, Board<W>)> {
        let m = x.erasable(j);
        let mut out = Vec::new();
        let mut sub = m;
        loop {
            for sign in [PLUS, MINUS] {
                out.push((sub, sign, x.apply(j, sub, sign, V::BURN)));
            }
            if V::FORCED || sub.is_zero() {
                break;
            }
            sub = sub.wrapping_sub(W::ONE) & m;
        }
        out
    }

    #[inline(always)]
    fn probe(&self, x: Board<W>, r: u32) -> Option<u32> {
        let h = hash_key(x.bits, x.len, r);
        self.table
            .find(h, |e| {
                let (bits, len, er) = (e.bits, e.len, e.r);
                bits == x.bits && len as u32 == x.len && er as u32 == r
            })
            .map(|e| e.val as u32)
    }

    #[inline(always)]
    fn store(&mut self, x: Board<W>, r: u32, val: u32) {
        if self.table.len() >= self.cfg.max_states {
            self.aborted = true;
            return;
        }
        let h = hash_key(x.bits, x.len, r);
        let e = Entry { bits: x.bits, len: x.len as u8, r: r as u8, val: val as u8 };
        self.table.insert_unique(h, e, |e| {
            let (bits, len, er) = (e.bits, e.len, e.r);
            hash_key(bits, len as u32, er as u32)
        });
    }

    /// `(bank, key)` with `V(x, r) = bank + V(key, r)`.
    #[inline(always)]
    pub fn reduce(&self, x: Board<W>) -> (u32, Board<W>) {
        if V::BURN && self.cfg.qv_reduce { x.reduce_qv() } else { x.reduce() }
    }

    #[inline(always)]
    fn child(&mut self, x: Board<W>, r: u32) -> u32 {
        let (bank, core) = self.reduce(x);
        bank + self.search(core, r)
    }

    /// `V(x, r)` for a banked, canonical core `x`.
    fn search(&mut self, x: Board<W>, r: u32) -> u32 {
        if self.aborted {
            return 0;
        }
        self.nodes += 1;
        let c = x.count();
        let empt = x.empties();
        if r == 0 || empt.is_zero() {
            return c;
        }
        if let Some(v) = self.probe(x, r) {
            return v;
        }
        let mut ub = c + r.min(empt.count_ones());
        let mut best = c;
        if self.cfg.monotone {
            if r >= 2
                && let Some(lb) = self.probe(x, r - 1)
            {
                best = best.max(lb);
            }
            if let Some(u) = self.probe(x, r + 1) {
                ub = ub.min(u);
            }
        }
        if best < ub {
            // Centre-out order: middle cells first.
            let mut order = [0u32; 64];
            let mut k = 0;
            for j in cells_of(empt) {
                order[k] = j;
                k += 1;
            }
            let mid2 = x.len as i32 - 1;
            order[..k].sort_unstable_by_key(|&j| (2 * j as i32 - mid2).abs());
            for &j in &order[..k] {
                let m = x.erasable(j);
                let lv = self.labeler(x, j, m, r, best);
                if lv > best {
                    best = lv;
                    if best >= ub {
                        break;
                    }
                }
            }
        }
        if !self.aborted {
            self.store(x, r, best);
        }
        best
    }

    /// An upper bound on `V(y, r)` that needs no search: the trivial bound
    /// `|y| + min(r, #empty)`, tightened by a stored value at `r` or, with
    /// monotone bounds on, at `r + 1`.
    #[inline(always)]
    fn cheap_upper(&self, y: Board<W>, r: u32) -> u32 {
        let mut ub = y.count() + r.min(y.empties().count_ones());
        let (bank, core) = self.reduce(y);
        if let Some(v) = self.probe(core, r) {
            return bank + v;
        }
        if self.cfg.monotone
            && let Some(v) = self.probe(core, r + 1)
        {
            ub = ub.min(bank + v);
        }
        ub
    }

    /// Minimum over the Labeler's options for pointer `j`, stopping early once a
    /// child is `<= alpha`.
    ///
    /// A first pass looks for a child whose cheap upper bound is already
    /// `<= alpha`. Such a child proves the minimum is `<= alpha` without any search,
    /// and the bound itself is returned (the caller only tests `> alpha`).
    #[inline(always)]
    fn labeler(&mut self, x: Board<W>, j: u32, m: W, r: u32, alpha: u32) -> u32 {
        let mut sub = m;
        loop {
            for sign in [PLUS, MINUS] {
                let ub = self.cheap_upper(x.apply(j, sub, sign, V::BURN), r - 1);
                if ub <= alpha {
                    return ub;
                }
            }
            if V::FORCED || sub.is_zero() {
                break;
            }
            sub = sub.wrapping_sub(W::ONE) & m;
        }

        let mut worst = u32::MAX;
        let mut sub = m;
        loop {
            for sign in [PLUS, MINUS] {
                let v = self.child(x.apply(j, sub, sign, V::BURN), r - 1);
                if v < worst {
                    worst = v;
                    if worst <= alpha {
                        return worst;
                    }
                }
            }
            // Forced variants have the single option E = m. The full variant walks
            // every submask of m, from m itself down to and including 0.
            if V::FORCED || sub.is_zero() {
                return worst;
            }
            sub = sub.wrapping_sub(W::ONE) & m;
        }
    }
}
