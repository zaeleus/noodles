//! Variant record info field value.

pub mod array;

use std::io;

pub use self::array::Array;

/// A variant record info field value.
#[derive(Debug)]
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

impl<'a> TryFrom<Value<'a>> for crate::variant::record_buf::info::field::Value {
    type Error = io::Error;

    fn try_from(value: Value<'a>) -> Result<Self, Self::Error> {
        match value {
            Value::Integer(n) => Ok(Self::Integer(n)),
            Value::Float(n) => Ok(Self::Float(n)),
            Value::Flag => Ok(Self::Flag),
            Value::Character(c) => Ok(Self::Character(c)),
            Value::String(s) => Ok(Self::String(s.into())),
            Value::Array(array) => array.try_into().map(Self::Array),
        }
    }
}
