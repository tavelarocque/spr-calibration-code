//! Exact solvers for SPR(n, s): the optimised memoised minimax ([`Solver`]), the
//! unreduced brute force used to validate it ([`brute::Brute`]), strategy
//! extraction and the erasure-lemma search.

pub mod brute;
pub mod erasure;
pub mod gap;
pub mod policy;
pub mod solver;
pub mod strategy;

pub use solver::{Config, Solver};
