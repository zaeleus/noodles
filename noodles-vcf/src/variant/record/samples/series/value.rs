//! Variant record samples value.

pub mod array;

pub use self::array::Array;

/// A variant record samples value.
#[derive(Debug)]
pub enum Value<'a> {
    /// A 32-bit integer.
    Integer(i32),
    /// A single-precision floating-point.
    Float(f32),
    /// A character.
    Character(char),
    /// A string.
    String(&'a str),
    /// An array.
    Array(Array<'a>),
}
