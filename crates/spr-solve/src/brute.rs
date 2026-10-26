//! Brute-force reference: memoised minimax on raw `Vec<u8>` boards, straight from
//! the definition. No canonicalisation, no banking, no cutoffs, no bounds. The
//! board model is `spr_core::naive`, which shares no code with the packed board.

use spr_core::VariantKind;
use spr_core::naive;
use std::collections::HashMap;

pub struct Brute {
    kind: VariantKind,
    memo: HashMap<(Vec<u8>, u32), u32>,
}

impl Brute {
    pub fn new(kind: VariantKind) -> Self {
        Brute { kind, memo: HashMap::new() }
    }

    /// `V(x, 0) = |x|`, `V(x, r) = max(|x|, max_j min_{E, sign} V(x', r - 1))`.
    pub fn value(&mut self, x: &[u8], r: u32) -> u32 {
        if r == 0 {
            return naive::count(x);
        }
        if let Some(&v) = self.memo.get(&(x.to_vec(), r)) {
            return v;
        }
        let mut best = naive::count(x);
        for j in naive::empties(x) {
            let mut worst = u32::MAX;
            for e in naive::erasure_options(&naive::erasable(x, j), self.kind) {
                for sign in [spr_core::PLUS, spr_core::MINUS] {
                    let y = naive::apply(x, j, &e, sign, self.kind);
                    worst = worst.min(self.value(&y, r - 1));
                }
            }
            best = best.max(worst);
        }
        self.memo.insert((x.to_vec(), r), best);
        best
    }

    pub fn opt(&mut self, n: usize, s: u32) -> u32 {
        self.value(&vec![0u8; n], s)
    }
}
