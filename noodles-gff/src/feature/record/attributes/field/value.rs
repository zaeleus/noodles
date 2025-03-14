//! Feature record attributes field value.

mod array;

use std::iter;

pub use self::array::Array;

/// A feature record attributes field value.
pub enum Value<'r> {
    /// A string.
    String(&'r str),
    /// An array.
    Array(Box<dyn Array + 'r>),
}

impl Value<'_> {
    /// Returns the value as a string, if the value is a string.
    pub fn as_string(&self) -> Option<&str> {
        match self {
            Self::String(value) => Some(value),
            Self::Array(_) => None,
        }
    }

    /// Returns the value as an array, if the value is an array.
    pub fn as_array(&self) -> Option<&dyn Array> {
        match self {
            Self::String(_) => None,
            Self::Array(array) => Some(array.as_ref()),
        }
    }

    /// Returns an iterator over values.
    pub fn iter(&self) -> Box<dyn Iterator<Item = &str> + '_> {
        match self {
            Self::String(value) => Box::new(iter::once(*value)),
            Self::Array(array) => Box::new(array.iter()),
        }
    }
}
