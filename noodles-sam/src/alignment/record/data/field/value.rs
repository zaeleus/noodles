//! Alignment record data field value.

pub mod array;

pub use self::array::Array;

/// An alignment record data field value.
pub enum Value<'a> {
    /// A character (`A`).
    Character(u8),
    /// An 8-bit integer (`c`).
    Int8(i8),
    /// An 8-bit unsigned integer (`C`).
    UInt8(u8),
    /// A 16-bit integer (`s`).
    Int16(i16),
    /// A 16-bit unsigned integer (`S`).
    UInt16(u16),
    /// A 32-bit integer (`i`).
    Int32(i32),
    /// A 32-bit unsigned integer (`I`).
    UInt32(u32),
    /// A single-precision floating-point (`f`).
    Float(f32),
    /// A string (`Z`).
    String(&'a [u8]),
    /// A hex string (`H`).
    Hex(&'a [u8]),
    /// An array (`B`).
    Array(Array<'a>),
}
