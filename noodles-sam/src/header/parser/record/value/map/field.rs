pub mod tag;

use std::{error, fmt, str};

pub use self::tag::parse_tag;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    MissingDelimiter,
    InvalidDelimiter,
    MissingSeparator,
    InvalidSeparator,
    InvalidValue(str::Utf8Error),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidValue(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingDelimiter => write!(f, "missing delimiter"),
            Self::InvalidDelimiter => write!(f, "invalid delimiter"),
            Self::MissingSeparator => write!(f, "missing separator"),
            Self::InvalidSeparator => write!(f, "invalid separator"),
            Self::InvalidValue(_) => write!(f, "invalid value"),
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

pub(super) fn parse_value<'a>(src: &mut &'a [u8]) -> Result<&'a str, ParseError> {
    use memchr::memchr;

    const DELIMITER: u8 = b'\t';

    let i = memchr(DELIMITER, src).unwrap_or(src.len());
    let (buf, rest) = src.split_at(i);

    *src = rest;

    str::from_utf8(buf).map_err(ParseError::InvalidValue)
}
