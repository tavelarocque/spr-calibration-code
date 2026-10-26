//! The oblivious no-reuse Pointer (DDFGKO Definition 5.1): float LP search for eps*,
//! exact rational certification, constructions (binomial cube, oblivious binary
//! search, products), and ov*. Definitions follow Dagan et al., arXiv:2406.13668, Definition 6.1.

pub mod constructions;
pub mod exact;
pub mod lp;
pub mod rational;
pub mod search;
pub mod strategy;
pub mod tuples;
