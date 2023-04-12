//! SAM record data field and components.

pub mod tag;
pub mod ty;
pub mod value;

pub use self::{tag::Tag, ty::Type, value::Value};

use std::{error, fmt};

/// An error returned when a raw SAM record data field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
    /// The data field tag is invalid.
    InvalidTag(tag::ParseError),
    /// The data field type is invalid.
    InvalidType(ty::ParseError),
    /// The data field value is invalid.
    InvalidValue,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidTag(e) => Some(e),
            Self::InvalidType(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidTag(_) => f.write_str("invalid tag"),
            Self::InvalidType(_) => f.write_str("invalid type"),
            Self::InvalidValue => f.write_str("invalid value"),
        }
    }
}
