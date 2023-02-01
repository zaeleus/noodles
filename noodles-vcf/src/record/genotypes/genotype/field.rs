//! VCF record genotype field.

pub mod value;

pub use self::value::Value;

use std::{error, fmt};

use crate::header::format::Key;

/// A VCF record genotype field.
#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    key: Key,
    value: Option<Value>,
}

/// An error returned when a raw VCF record info field value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The value is invalid.
    InvalidValue(value::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidValue(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidValue(_) => f.write_str("invalid value"),
        }
    }
}
