pub mod tag;
pub mod value;

use std::{error, fmt};

pub use self::{tag::parse_tag, value::parse_value};

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

    if let Some((b, rest)) = src.split_first() {
        if *b == DELIMITER {
            *src = rest;
            Ok(())
        } else {
            Err(ParseError::InvalidDelimiter)
        }
    } else {
        Err(ParseError::MissingDelimiter)
    }
}

pub(super) fn consume_separator(src: &mut &[u8]) -> Result<(), ParseError> {
    const SEPARATOR: u8 = b':';

    if let Some((b, rest)) = src.split_first() {
        if *b == SEPARATOR {
            *src = rest;
            Ok(())
        } else {
            Err(ParseError::InvalidSeparator)
        }
    } else {
        Err(ParseError::MissingSeparator)
    }
}
