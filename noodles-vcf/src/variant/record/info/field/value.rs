//! Variant record info field value.

pub mod array;

pub use self::array::Array;

/// A variant record info field value.
pub enum Value<'a> {
    /// A 32-bit integer.
    Integer(i32),
    /// A single-precision floating-point.
    Float(f32),
    /// A boolean.
    Flag,
    /// A character.
    Character(char),
    /// A string.
    String(&'a str),
    /// An array.
    Array(Array<'a>),
}
