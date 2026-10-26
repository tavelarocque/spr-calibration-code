//! Unreduced reference implementation on `Vec<u8>` boards, written to mirror
//! reference/spr_ref.py line by line. Used only to test the packed code and as
//! the board model of the brute-force solver. No packing, no canon, no banking.

use crate::board::{BURNED, MINUS, PLUS};
use crate::variant::VariantKind;

pub fn count(x: &[u8]) -> u32 {
    x.iter().filter(|&&c| c == PLUS || c == MINUS).count() as u32
}

pub fn empties(x: &[u8]) -> Vec<usize> {
    (0..x.len()).filter(|&j| x[j] == 0).collect()
}

pub fn erasable(x: &[u8], j: usize) -> Vec<usize> {
    (0..x.len()).filter(|&i| (i < j && x[i] == MINUS) || (i > j && x[i] == PLUS)).collect()
}

/// Reverse the strip and swap '+' with '-'.
pub fn reflect(x: &[u8]) -> Vec<u8> {
    x.iter()
        .rev()
        .map(|&c| match c {
            PLUS => MINUS,
            MINUS => PLUS,
            c => c,
        })
        .collect()
}

/// Every erasure set the Labeler may choose, as lists of cells.
pub fn erasure_options(er: &[usize], kind: VariantKind) -> Vec<Vec<usize>> {
    if kind.forced() {
        return vec![er.to_vec()];
    }
    (0u64..1 << er.len())
        .map(|m| (0..er.len()).filter(|&b| m >> b & 1 == 1).map(|b| er[b]).collect())
        .collect()
}

/// Successor board: erase `e`, write `sign` at `j`.
pub fn apply(x: &[u8], j: usize, e: &[usize], sign: u8, kind: VariantKind) -> Vec<u8> {
    let mut y = x.to_vec();
    for &i in e {
        y[i] = if kind.burn() { BURNED } else { 0 };
    }
    y[j] = sign;
    y
}

/// Every successor `(j, E, sign, x')` of `x`.
pub fn successors(x: &[u8], kind: VariantKind) -> Vec<(usize, Vec<usize>, u8, Vec<u8>)> {
    let mut out = Vec::new();
    for j in empties(x) {
        for e in erasure_options(&erasable(x, j), kind) {
            for sign in [PLUS, MINUS] {
                let y = apply(x, j, &e, sign, kind);
                out.push((j, e.clone(), sign, y));
            }
        }
    }
    out
}

/// Length of the maximal '+' prefix and '-' suffix, by direct scan.
pub fn bank_lengths(x: &[u8]) -> (usize, usize) {
    let p = x.iter().take_while(|&&c| c == PLUS).count();
    let q = x[p..].iter().rev().take_while(|&&c| c == MINUS).count();
    (p, q)
}

/// QV reduction by direct scan: drop burned cells, then drop (and count) every
/// '+' left of the first empty cell and every '-' right of the last empty cell.
/// With no empty cell, every sign is dropped and counted.
pub fn qv_reduce(x: &[u8]) -> (u32, Vec<u8>) {
    let y: Vec<u8> = x.iter().copied().filter(|&c| c != BURNED).collect();
    let Some(first) = y.iter().position(|&c| c == 0) else {
        return (count(&y), Vec::new());
    };
    let last = y.iter().rposition(|&c| c == 0).unwrap();
    let mut bank = 0;
    let mut core = Vec::new();
    for (i, &c) in y.iter().enumerate() {
        if (c == PLUS && i < first) || (c == MINUS && i > last) {
            bank += 1;
        } else {
            core.push(c);
        }
    }
    (bank, core)
}
