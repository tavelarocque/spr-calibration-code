//! Ordered s-tuples of distinct cells in [n] (cells are 0-based here), their prefixes
//! and the reflection `c -> n - 1 - c`.
//!
//! Prefixes of length `1 ..= s - 1` are ranked densely: a prefix of length `i` gets
//! `offset[i] + rank`, where `rank` is its mixed-radix rank among the `n^(i)` injective
//! sequences of length `i` (falling factorial). The rank of `(c_1, ..., c_i)` uses, at
//! position `j`, the number of still-unused cells below `c_j`.

/// Falling factorial `n (n - 1) ... (n - k + 1)`.
pub fn falling(n: usize, k: usize) -> usize {
    (0..k).map(|j| n - j).product()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Space {
    pub n: usize,
    pub s: usize,
}

impl Space {
    pub fn new(n: usize, s: usize) -> Self {
        assert!(n <= 128, "cells are tracked in a u128 mask");
        assert!(s >= 1 && s <= n, "need 1 <= s <= n");
        Space { n, s }
    }

    pub fn num_tuples(&self) -> usize {
        falling(self.n, self.s)
    }

    /// Number of prefixes of length 1 ..= s - 1 (the constraint rows).
    pub fn num_prefixes(&self) -> usize {
        (1..self.s).map(|i| falling(self.n, i)).sum()
    }

    /// Row offset of prefixes of length `i` (1-based length).
    pub fn offset(&self, i: usize) -> usize {
        (1..i).map(|l| falling(self.n, l)).sum()
    }

    /// Dense index of the prefix `p` (length 1 ..= s - 1).
    pub fn prefix_index(&self, p: &[u8]) -> usize {
        let i = p.len();
        debug_assert!(i >= 1 && i < self.s);
        let mut used: u128 = 0;
        let mut rank = 0usize;
        for (j, &c) in p.iter().enumerate() {
            let below = (used & ((1u128 << c) - 1)).count_ones() as usize;
            let idx = c as usize - below;
            // Remaining positions after j within the prefix: i - j - 1, each with
            // n - j - 1, n - j - 2, ... choices.
            rank = rank * (self.n - j) + idx;
            used |= 1u128 << c;
        }
        self.offset(i) + rank
    }

    /// Inverse of [`Space::prefix_index`].
    pub fn prefix_from_index(&self, mut idx: usize) -> Vec<u8> {
        let mut i = 1;
        while idx >= falling(self.n, i) {
            idx -= falling(self.n, i);
            i += 1;
        }
        let mut digits = vec![0usize; i];
        for j in (0..i).rev() {
            let radix = self.n - j;
            digits[j] = idx % radix;
            idx /= radix;
        }
        let mut used: u128 = 0;
        let mut out = Vec::with_capacity(i);
        for &d in &digits {
            // The d-th unused cell.
            let mut count = 0;
            let mut c = 0u8;
            loop {
                if used >> c & 1 == 0 {
                    if count == d {
                        break;
                    }
                    count += 1;
                }
                c += 1;
            }
            used |= 1u128 << c;
            out.push(c);
        }
        out
    }

    pub fn reflect(&self, t: &[u8]) -> Vec<u8> {
        t.iter().map(|&c| (self.n - 1) as u8 - c).collect()
    }

    /// All tuples, lexicographically, flattened with stride `s`.
    pub fn all_tuples(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.num_tuples() * self.s);
        let mut cur = Vec::with_capacity(self.s);
        fn rec(n: usize, s: usize, used: u128, cur: &mut Vec<u8>, out: &mut Vec<u8>) {
            if cur.len() == s {
                out.extend_from_slice(cur);
                return;
            }
            for c in 0..n as u8 {
                if used >> c & 1 == 0 {
                    cur.push(c);
                    rec(n, s, used | 1u128 << c, cur, out);
                    cur.pop();
                }
            }
        }
        rec(self.n, self.s, 0, &mut cur, &mut out);
        out
    }
}

/// For a tuple `t`, `(up[i], down[i])` for each position `i` (0-based): every later
/// cell is strictly greater (up) or strictly smaller (down) than `t[i]`. For the last
/// position both are vacuously true.
pub fn up_down(t: &[u8]) -> (Vec<bool>, Vec<bool>) {
    let s = t.len();
    let mut up = vec![true; s];
    let mut down = vec![true; s];
    let (mut mn, mut mx) = (u8::MAX, 0u8);
    let mut any = false;
    for i in (0..s).rev() {
        if any {
            up[i] = mn > t[i];
            down[i] = mx < t[i];
        }
        mn = mn.min(t[i]);
        mx = mx.max(t[i]);
        any = true;
    }
    (up, down)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prefix_ranks_are_a_bijection() {
        for (n, s) in [(5, 3), (7, 4), (6, 6)] {
            let sp = Space::new(n, s);
            let mut seen = vec![false; sp.num_prefixes()];
            let all = sp.all_tuples();
            for t in all.chunks(s) {
                for i in 1..s {
                    let idx = sp.prefix_index(&t[..i]);
                    assert_eq!(sp.prefix_from_index(idx), t[..i].to_vec());
                    seen[idx] = true;
                }
            }
            assert!(seen.iter().all(|&b| b));
            assert_eq!(all.len(), sp.num_tuples() * s);
        }
    }

    #[test]
    fn up_down_matches_definition() {
        let t = [3u8, 5, 1, 4];
        let (u, d) = up_down(&t);
        assert_eq!(u, vec![false, false, true, true]);
        assert_eq!(d, vec![false, true, false, true]);
    }
}
