//! Variant record info field array value.

mod values;

pub use self::values::Values;

/// A variant record info field array value.
pub enum Array<'a> {
    /// A 32-bit integer array.
    Integer(Box<dyn Values<'a, i32> + 'a>),
    /// A single-precision floating-point array..
    Float(Box<dyn Values<'a, f32> + 'a>),
    /// A character array.
    Character(Box<dyn Values<'a, char> + 'a>),
    /// A string array.
    String(Box<dyn Values<'a, &'a str> + 'a>),
}
