//! Machine words that hold a packed board, 2 bits per cell.
//!
//! [`Word`] is implemented for `u64` (up to 32 cells) and `u128` (up to 64 cells).
//! Everything the board code needs is expressed through this trait so that the
//! solver is monomorphised per width and pays nothing for the abstraction.

use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{BitAnd, BitOr, BitXor, Not, Shl, Shr};

pub trait Word:
    Copy
    + Eq
    + Ord
    + Hash
    + Debug
    + Default
    + BitAnd<Output = Self>
    + BitOr<Output = Self>
    + BitXor<Output = Self>
    + Not<Output = Self>
    + Shl<u32, Output = Self>
    + Shr<u32, Output = Self>
    + Send
    + Sync
    + 'static
{
    /// Width in bits.
    const BITS: u32;
    /// Maximum number of cells, `BITS / 2`.
    const CELLS: u32 = Self::BITS / 2;
    const ZERO: Self;
    const ONE: Self;
    const ONES: Self;
    /// Bit 0 of every cell: `0b0101...01`.
    const EVEN: Self;

    fn count_ones(self) -> u32;
    fn trailing_zeros(self) -> u32;
    fn leading_zeros(self) -> u32;
    fn reverse_bits(self) -> Self;
    fn wrapping_sub(self, rhs: Self) -> Self;
    /// The two 64-bit halves, low first, for hashing.
    fn halves(self) -> (u64, u64);
    fn from_u64(v: u64) -> Self;

    /// Mask of the low `bits` bits. Valid for `bits` in `0..=BITS`.
    #[inline(always)]
    fn low_bits(bits: u32) -> Self {
        if bits >= Self::BITS {
            Self::ONES
        } else {
            (Self::ONE << bits).wrapping_sub(Self::ONE)
        }
    }

    #[inline(always)]
    fn is_zero(self) -> bool {
        self == Self::ZERO
    }
}

macro_rules! impl_word {
    ($t:ty, $even:expr) => {
        impl Word for $t {
            const BITS: u32 = <$t>::BITS;
            const ZERO: Self = 0;
            const ONE: Self = 1;
            const ONES: Self = !0;
            const EVEN: Self = $even;

            #[inline(always)]
            fn count_ones(self) -> u32 {
                <$t>::count_ones(self)
            }
            #[inline(always)]
            fn trailing_zeros(self) -> u32 {
                <$t>::trailing_zeros(self)
            }
            #[inline(always)]
            fn leading_zeros(self) -> u32 {
                <$t>::leading_zeros(self)
            }
            #[inline(always)]
            fn reverse_bits(self) -> Self {
                <$t>::reverse_bits(self)
            }
            #[inline(always)]
            fn wrapping_sub(self, rhs: Self) -> Self {
                <$t>::wrapping_sub(self, rhs)
            }
            #[inline(always)]
            fn halves(self) -> (u64, u64) {
                (self as u64, ((self as u128) >> 64) as u64)
            }
            #[inline(always)]
            fn from_u64(v: u64) -> Self {
                v as $t
            }
        }
    };
}

impl_word!(u64, 0x5555_5555_5555_5555);
impl_word!(u128, 0x5555_5555_5555_5555_5555_5555_5555_5555);
