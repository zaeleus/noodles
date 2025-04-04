//! Feature record attributes field value.

mod array;

use std::{borrow::Cow, io, iter};

use bstr::BStr;

pub use self::array::Array;

/// A feature record attributes field value.
pub enum Value<'r> {
    /// A string.
    String(Cow<'r, BStr>),
    /// An array.
    Array(Box<dyn Array<'r> + 'r>),
}

impl<'r> Value<'r> {
    /// Returns the value as a string, if the value is a string.
    pub fn as_string(&self) -> Option<&BStr> {
        match self {
            Self::String(value) => Some(value.as_ref()),
            Self::Array(_) => None,
        }
    }

    /// Returns the value as an array, if the value is an array.
    pub fn as_array(&self) -> Option<&'r dyn Array> {
        match self {
            Self::String(_) => None,
            Self::Array(array) => Some(array.as_ref()),
        }
    }

    /// Returns an iterator over values.
    pub fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Cow<'r, BStr>>> + '_> {
        match self {
            Self::String(value) => Box::new(iter::once(Ok(value.clone()))),
            Self::Array(array) => Box::new(array.iter()),
        }
    }
}
