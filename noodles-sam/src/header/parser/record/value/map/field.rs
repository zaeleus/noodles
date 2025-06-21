pub mod tag;
pub mod value;

use std::{error, fmt};

pub use self::{tag::parse_tag, value::parse_value};
use crate::header::parser::record::split_off_first;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    MissingDelimiter,
    InvalidDelimiter,
    MissingSeparator,
    InvalidSeparator,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingDelimiter => write!(f, "missing delimiter"),
            Self::InvalidDelimiter => write!(f, "invalid delimiter"),
            Self::MissingSeparator => write!(f, "missing separator"),
            Self::InvalidSeparator => write!(f, "invalid separator"),
        }
    }
}

pub(super) fn consume_delimiter(src: &mut &[u8]) -> Result<(), ParseError> {
    const DELIMITER: u8 = b'\t';

    match split_off_first(src) {
        Some(&DELIMITER) => Ok(()),
        Some(_) => Err(ParseError::InvalidDelimiter),
        None => Err(ParseError::MissingDelimiter),
    }
}

pub(super) fn consume_separator(src: &mut &[u8]) -> Result<(), ParseError> {
    const SEPARATOR: u8 = b':';

    match split_off_first(src) {
        Some(&SEPARATOR) => Ok(()),
        Some(_) => Err(ParseError::InvalidSeparator),
        None => Err(ParseError::MissingSeparator),
    }
}
