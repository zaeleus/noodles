//! SAM record data field value subtype.

/// A SAM record data field value subtype.
///
/// Only arrays have subtypes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Subtype {
    /// 8-bit integer (`c`).
    Int8,
    /// 8-bit unsigned integer (`C`).
    UInt8,
    /// 16-bit integer (`s`).
    Int16,
    /// 16-bit unsigned integer (`S`).
    UInt16,
    /// 32-bit integer (`i`).
    Int32,
    /// 32-bit unsigned integer (`I`).
    UInt32,
    /// Single-precision floating-point (`f`).
    Float,
}
