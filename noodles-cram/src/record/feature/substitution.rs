//! CRAM record substitution feature.

mod base;

pub use self::base::Base;

/// A substitution feature value.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Value {
    /// A substitution code used when reading and writing.
    Code(u8),
    /// A substitution-read base pair used when staging the record for writing.
    Bases(Base, Base),
}
