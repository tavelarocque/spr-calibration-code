//! The three game variants.
//!
//! A trait with associated constants rather than a runtime enum: the solver is
//! generic over `V: Variant`, so each variant gets its own compiled search loop
//! and the per-node checks below fold away at compile time. The CLI holds the
//! only runtime choice ([`VariantKind`]) and dispatches once, at the top.

/// Rules shared by every variant: the Pointer points at an empty (code 00) cell,
/// the Labeler erases a subset of the erasable set and writes '+' or '-'.
pub trait Variant: Copy + Default + Send + Sync + 'static {
    const KIND: VariantKind;
    /// Labeler must erase the whole erasable set.
    const FORCED: bool;
    /// Erased cells become burned (code 11) instead of empty.
    const BURN: bool;
}

#[derive(Clone, Copy, Default, Debug)]
pub struct Full;
#[derive(Clone, Copy, Default, Debug)]
pub struct Forced;
#[derive(Clone, Copy, Default, Debug)]
pub struct Qv;

impl Variant for Full {
    const KIND: VariantKind = VariantKind::Full;
    const FORCED: bool = false;
    const BURN: bool = false;
}

impl Variant for Forced {
    const KIND: VariantKind = VariantKind::Forced;
    const FORCED: bool = true;
    const BURN: bool = false;
}

impl Variant for Qv {
    const KIND: VariantKind = VariantKind::Qv;
    const FORCED: bool = true;
    const BURN: bool = true;
}

#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub enum VariantKind {
    Full,
    Forced,
    Qv,
}

impl VariantKind {
    pub const ALL: [VariantKind; 3] = [VariantKind::Full, VariantKind::Forced, VariantKind::Qv];

    pub fn name(self) -> &'static str {
        match self {
            VariantKind::Full => "full",
            VariantKind::Forced => "forced",
            VariantKind::Qv => "qv",
        }
    }

    pub fn forced(self) -> bool {
        !matches!(self, VariantKind::Full)
    }

    pub fn burn(self) -> bool {
        matches!(self, VariantKind::Qv)
    }

    /// The oracle's `make_solver(game, lmode)` arguments.
    pub fn oracle_args(self) -> (&'static str, &'static str) {
        match self {
            VariantKind::Full => ("spr", "any"),
            VariantKind::Forced => ("spr", "all"),
            VariantKind::Qv => ("qv", "all"),
        }
    }
}

impl std::str::FromStr for VariantKind {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "full" => Ok(VariantKind::Full),
            "forced" => Ok(VariantKind::Forced),
            "qv" => Ok(VariantKind::Qv),
            _ => Err(format!("unknown variant {s:?} (expected full, forced or qv)")),
        }
    }
}

impl std::fmt::Display for VariantKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}
