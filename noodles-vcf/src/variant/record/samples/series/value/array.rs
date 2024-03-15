//! Variant record samples array value.

mod values;

use std::{fmt, io};

pub use self::values::Values;

/// A variant record samples array value.
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

impl<'a> fmt::Debug for Array<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Integer(values) => f.debug_list().entries(values.iter()).finish(),
            Self::Float(values) => f.debug_list().entries(values.iter()).finish(),
            Self::Character(values) => f.debug_list().entries(values.iter()).finish(),
            Self::String(values) => f.debug_list().entries(values.iter()).finish(),
        }
    }
}

impl<'a> TryFrom<Array<'a>> for crate::variant::record_buf::samples::sample::value::Array {
    type Error = io::Error;

    fn try_from(array: Array<'a>) -> Result<Self, Self::Error> {
        match array {
            Array::Integer(values) => values.iter().collect::<io::Result<_>>().map(Self::Integer),
            Array::Float(values) => values.iter().collect::<io::Result<_>>().map(Self::Float),
            Array::Character(values) => values
                .iter()
                .collect::<io::Result<_>>()
                .map(Self::Character),
            Array::String(values) => values
                .iter()
                .map(|result| result.map(|value| value.map(String::from)))
                .collect::<io::Result<_>>()
                .map(Self::String),
        }
    }
}
