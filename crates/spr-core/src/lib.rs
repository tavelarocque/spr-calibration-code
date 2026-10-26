//! Game definition for Sign Preservation with Reuse: packed boards, move
//! generation, the three variants, and the two exact reductions (reflection and
//! banking) with their proofs in [`board`].

pub mod board;
pub mod naive;
pub mod variant;
pub mod word;

pub use board::{BURNED, Board, EMPTY, MINUS, PLUS, cell_char, cells_of};
pub use variant::{Forced, Full, Qv, Variant, VariantKind};
pub use word::Word;
