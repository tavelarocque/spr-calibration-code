//! Packed board: 2 bits per cell, cell `i` at bits `2i` and `2i + 1`.
//!
//! Codes: `00` empty, `01` '+', `10` '-', `11` burned (QV only: erased, never
//! pointed at again). These are the oracle's codes 0, 1, 2, 3. Bits above cell
//! `len - 1` are always zero.
//!
//! Most queries return a *spread mask*: a word with bit `2i` set iff cell `i` has
//! some property, all odd bits clear. Spread masks need no bit gathering, and a
//! spread mask `m` widens to "both bits of those cells" as `m | m << 1`.
//!
//! # Reflection (proof of correctness)
//!
//! Let `rho(i) = n - 1 - i` and let `sigma` swap '+' and '-', fixing empty and
//! burned. Define the reflection `y_i = sigma(x_{rho(i)})`. Then `V(x, r) = V(y, r)`
//! in every variant.
//!
//! Proof by induction on `r`. `|y| = |x|`, which settles `r = 0` and the stop
//! option. Cell `j` is empty in `x` iff `rho(j)` is empty in `y`, and burned cells
//! map to burned cells, so the Pointer's options correspond. For the erasable sets,
//! `rho` maps `{ i < j : x_i = - }` onto `{ i' > rho(j) : y_i' = + }` and
//! `{ i > j : x_i = + }` onto `{ i' < rho(j) : y_i' = - }`, so `rho(E_x(j)) =
//! E_y(rho(j))`. This preserves "any subset" and "the whole set". Erasing `E` and
//! writing `s` at `j` in `x` therefore corresponds to erasing `rho(E)` and writing
//! `sigma(s)` at `rho(j)` in `y`. The two successors are again reflections of each
//! other, so their values agree at `r - 1`, and the max/min over corresponding
//! options agree at `r`.
//!
//! Bit level: `reverse_bits` sends bit `2i + b` (b in {0, 1}) to
//! `BITS - 1 - 2i - b = 2(CELLS - 1 - i) + (1 - b)`. It moves cell `i` to cell
//! `CELLS - 1 - i` and swaps the two bits inside the cell, which maps codes
//! 00, 01, 10, 11 to 00, 10, 01, 11. That is exactly `sigma`. Shifting right by
//! `BITS - 2n` then moves cell `CELLS - 1 - i` to `n - 1 - i = rho(i)`. The zero
//! cells at positions `>= n` land below bit `BITS - 2n` and are shifted out.
//!
//! # Banking (proof of correctness)
//!
//! Suppose `x = P y S` where `P` is `p` cells all '+', `S` is `q` cells all '-'
//! and `y` is any board of length `m`. Then `V(x, r) = p + q + V(y, r)`, with the
//! same rules (reuse or not, any subset or forced) on both sides.
//!
//! Lemma: no cell of `P` is ever erased. Suppose one is, and take the first time
//! it happens, at cell `i` in `P` with the Pointer at `j`. A '+' at `i` is erasable
//! only if `j < i`, so `j` is in `P`. No cell of `P` has been erased yet, so `x_j`
//! is still '+' and `j` is not empty, a contradiction. Reuse makes no difference:
//! a cell can only become pointable again by being erased, which the lemma rules
//! out. Mirror argument for `S` (a '-' at `i` is erasable only if `j > i`). So the
//! Pointer only ever points into `y`, `E_x(j)` contains no cell of `P` or `S`,
//! and `E_x(j)`, shifted by `p`, is the erasable set of `y` at `j - p`.
//!
//! Induction on `r`: `|x| = p + q + |y|` settles `r = 0`, the stop option, and the
//! "no empty cell" case (`x` has an empty cell iff `y` does). For `r >= 1` the
//! Pointer's options in `x` and `y` correspond as above, and every successor of
//! `x` has the form `P y' S` with `y'` the corresponding successor of `y`. By
//! induction `V(P y' S, r - 1) = p + q + V(y', r - 1)`, and the constant `p + q`
//! passes through the min and max.
//!
//! [`Board::bank`] strips the *maximal* '+' prefix and '-' suffix, a special case.
//! In QV a burned cell is not '+', so it ends the prefix. That is sound (the lemma
//! only needs `P` to be all '+') but not maximal.
//!
//! Banking and reflection commute: reflection maps a '+' prefix of length `p` to a
//! '-' suffix of length `p` and vice versa, so `bank(reflect(x))` is the reflection
//! of `bank(x)` with the same bank. The memo key is therefore
//! `canon(bank(x))`, and `V(x, r) = bank + V(canon core, r)`.

use crate::word::Word;
use std::fmt;

pub const EMPTY: u8 = 0;
pub const PLUS: u8 = 1;
pub const MINUS: u8 = 2;
pub const BURNED: u8 = 3;

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub struct Board<W: Word> {
    pub bits: W,
    pub len: u32,
}

impl<W: Word> Board<W> {
    /// The all-empty board `0^n`.
    pub fn empty(n: u32) -> Self {
        assert!(n <= W::CELLS, "board of {n} cells does not fit in {} bits", W::BITS);
        Board { bits: W::ZERO, len: n }
    }

    /// Build from oracle cell codes (0, 1, 2, 3).
    pub fn from_cells(cells: &[u8]) -> Self {
        let mut b = Self::empty(cells.len() as u32);
        for (i, &c) in cells.iter().enumerate() {
            assert!(c <= BURNED);
            b.bits = b.bits | (W::from_u64(c as u64) << (2 * i as u32));
        }
        b
    }

    pub fn to_cells(&self) -> Vec<u8> {
        (0..self.len).map(|i| self.cell(i)).collect()
    }

    /// Parse `0`, `+`, `-`, `#` (burned).
    pub fn parse(s: &str) -> Option<Self> {
        let cells: Option<Vec<u8>> = s
            .chars()
            .map(|c| match c {
                '0' => Some(EMPTY),
                '+' => Some(PLUS),
                '-' => Some(MINUS),
                '#' => Some(BURNED),
                _ => None,
            })
            .collect();
        let cells = cells?;
        (cells.len() as u32 <= W::CELLS).then(|| Self::from_cells(&cells))
    }

    #[inline(always)]
    pub fn cell(&self, i: u32) -> u8 {
        let (lo, _) = ((self.bits >> (2 * i)) & W::from_u64(3)).halves();
        lo as u8
    }

    /// Spread mask of all cells `0..len`.
    #[inline(always)]
    pub fn valid(&self) -> W {
        W::EVEN & W::low_bits(2 * self.len)
    }

    #[inline(always)]
    fn lo(&self) -> W {
        self.bits & W::EVEN
    }

    #[inline(always)]
    fn hi(&self) -> W {
        (self.bits >> 1) & W::EVEN
    }

    /// Spread mask of '+' cells.
    #[inline(always)]
    pub fn plus(&self) -> W {
        self.lo() & !self.hi()
    }

    /// Spread mask of '-' cells.
    #[inline(always)]
    pub fn minus(&self) -> W {
        self.hi() & !self.lo()
    }

    /// Spread mask of empty (pointable) cells. Burned cells are not included.
    #[inline(always)]
    pub fn empties(&self) -> W {
        !(self.lo() | self.hi()) & self.valid()
    }

    /// Spread mask of burned cells.
    #[inline(always)]
    pub fn burned(&self) -> W {
        self.lo() & self.hi()
    }

    /// `|x|`, the number of signed cells.
    #[inline(always)]
    pub fn count(&self) -> u32 {
        (self.lo() ^ self.hi()).count_ones()
    }

    /// Spread mask of `{ i < j : x_i = - } union { i > j : x_i = + }`.
    #[inline(always)]
    pub fn erasable(&self, j: u32) -> W {
        let below = W::low_bits(2 * j);
        let above = !W::low_bits(2 * j + 2);
        (self.minus() & below) | (self.plus() & above)
    }

    /// Erase the cells of spread mask `erase` (to empty, or to burned if `burn`),
    /// then write `sign` into cell `j`, which must be empty.
    #[inline(always)]
    pub fn apply(&self, j: u32, erase: W, sign: u8, burn: bool) -> Self {
        debug_assert_eq!(self.cell(j), EMPTY);
        debug_assert!((erase & !self.erasable(j)).is_zero());
        let wide = erase | (erase << 1);
        let mut bits = self.bits & !wide;
        if burn {
            bits = bits | wide;
        }
        bits = bits | (W::from_u64(sign as u64) << (2 * j));
        Board { bits, len: self.len }
    }

    /// Reverse the strip and swap '+' with '-'. See the module docs for the proof.
    #[inline(always)]
    pub fn reflect(&self) -> Self {
        if self.len == 0 {
            return *self;
        }
        Board { bits: self.bits.reverse_bits() >> (W::BITS - 2 * self.len), len: self.len }
    }

    /// The numerically smaller of `x` and its reflection.
    #[inline(always)]
    pub fn canon(&self) -> Self {
        let r = self.reflect();
        if r.bits < self.bits { r } else { *self }
    }

    /// Strip the maximal '+' prefix and '-' suffix. Returns `(p + q, core)`.
    #[inline(always)]
    pub fn bank(&self) -> (u32, Self) {
        // Cells at index >= len are 00, so they stop the prefix scan by themselves.
        // When every one of the CELLS cells is '+', trailing_zeros returns BITS.
        let p = (W::EVEN & !self.plus()).trailing_zeros() / 2;
        let p = p.min(self.len);
        let not_minus = self.valid() & !self.minus();
        let q = if not_minus.is_zero() {
            self.len
        } else {
            let top_cell = (W::BITS - 1 - not_minus.leading_zeros()) / 2;
            self.len - 1 - top_cell
        };
        if p + q >= self.len {
            // Only possible when the board is '+'^p '-'^q exactly.
            return (self.len, Board { bits: W::ZERO, len: 0 });
        }
        let m = self.len - p - q;
        let bits = (self.bits >> (2 * p)) & W::low_bits(2 * m);
        (p + q, Board { bits, len: m })
    }

    /// Bank, then canonicalise the core. `V(x, r) = bank + V(core, r)`.
    #[inline(always)]
    pub fn reduce(&self) -> (u32, Self) {
        let (b, core) = self.bank();
        (b, core.canon())
    }

    /// Delete the cells of spread mask `m`, closing the gaps. Cells keep their order.
    pub fn delete_cells(&self, mut m: W) -> Self {
        let (mut bits, mut len) = (self.bits, self.len);
        // Highest first, so the indices of the remaining marked cells stay valid.
        while !m.is_zero() {
            let top = (W::BITS - 1 - m.leading_zeros()) / 2;
            let low = bits & W::low_bits(2 * top);
            let high = if 2 * top + 2 >= W::BITS { W::ZERO } else { (bits >> (2 * top + 2)) << (2 * top) };
            bits = low | high;
            len -= 1;
            m = m & !(W::ONE << (2 * top));
        }
        Board { bits, len }
    }

    /// QV only: delete burned cells and permanent signs, bank the permanent signs,
    /// canonicalise. `V_qv(x, r) = bank + V_qv(core, r)`.
    ///
    /// # Proof
    ///
    /// Deletion lemma: if a cell `c` is never pointed at and never erased in any
    /// continuation from `x`, then `V(x, r) = [c signed] + V(x \ c, r)`, where
    /// `x \ c` closes the gap. The Pointer's options are the empty cells other than
    /// `c`, in the same order. Erasable sets are defined by the position of each
    /// signed cell relative to `j`, and deletion preserves relative order, so
    /// `E_x(j)` is `E_{x \ c}(j)` relabelled (it never contains `c`, since QV erases
    /// the whole erasable set and `c` is never erased). Successors correspond
    /// (`x' \ c = (x \ c)'`), `|x| = [c signed] + |x \ c|`, and induction on `r`
    /// gives the identity, exactly as for banking.
    ///
    /// Burned cells satisfy the hypothesis: they are not empty, so never pointed
    /// at, and not signed, so never in an erasable set. They count 0.
    ///
    /// Permanent signs satisfy it too. QV never creates an empty cell (erased cells
    /// become burned), so the set of empty cells only shrinks. A '+' at `i` is
    /// erased only when the Pointer picks an empty `j < i`. If no cell left of `i`
    /// is empty now, none ever will be, so the '+' is never erased, and as a sign
    /// it is never pointed at. Mirror for a '-' with no empty cell to its right.
    /// Each counts 1. If the board has no empty cell at all, every sign is
    /// permanent and the core is the empty board.
    ///
    /// Deleting cells one at a time keeps the hypothesis for the others (deletion
    /// does not change which cells are empty, or which side of them a sign is on),
    /// so all of them can be removed together. Reflection maps the QV-reduced form
    /// of `x` to that of its reflection, so canonicalising afterwards is sound.
    pub fn reduce_qv(&self) -> (u32, Self) {
        let b = self.delete_cells(self.burned());
        let e = b.empties();
        if e.is_zero() {
            return (b.count(), Board { bits: W::ZERO, len: 0 });
        }
        let first = e.trailing_zeros() / 2;
        let last = (W::BITS - 1 - e.leading_zeros()) / 2;
        let perm = (b.plus() & W::low_bits(2 * first)) | (b.minus() & !W::low_bits(2 * last + 2));
        let core = b.delete_cells(perm);
        (perm.count_ones(), core.canon())
    }

    /// True iff every '+' lies strictly left of every '-'.
    pub fn plus_left_of_minus(&self) -> bool {
        let p = self.plus();
        let m = self.minus();
        if p.is_zero() || m.is_zero() {
            return true;
        }
        let last_plus = W::BITS - 1 - p.leading_zeros();
        let first_minus = m.trailing_zeros();
        last_plus < first_minus
    }

    /// Widen to another word size.
    pub fn convert<V: Word>(&self) -> Board<V> {
        Board::<V>::from_cells(&self.to_cells())
    }
}

/// Iterate the cell indices of a spread mask, lowest first.
#[inline(always)]
pub fn cells_of<W: Word>(mut m: W) -> impl Iterator<Item = u32> {
    std::iter::from_fn(move || {
        if m.is_zero() {
            None
        } else {
            let tz = m.trailing_zeros();
            m = m & m.wrapping_sub(W::ONE);
            Some(tz / 2)
        }
    })
}

pub fn cell_char(c: u8) -> char {
    match c {
        EMPTY => '0',
        PLUS => '+',
        MINUS => '-',
        _ => '#',
    }
}

impl<W: Word> fmt::Display for Board<W> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.len {
            write!(f, "{}", cell_char(self.cell(i)))?;
        }
        Ok(())
    }
}
