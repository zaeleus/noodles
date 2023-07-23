#![allow(dead_code)]

mod kind;
mod value;

use std::{error, fmt};

use self::{kind::parse_kind, value::parse_value};
use crate::header::Record;

/// An error returned when a SAM header record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The header prefix is missing.
    MissingPrefix,
    /// The input is invalid.
    InvalidKind(kind::ParseError),
    /// The delimiter is missing.
    MissingDelimiter,
    /// The value is invalid.
    InvalidValue(value::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKind(e) => Some(e),
            Self::InvalidValue(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingPrefix => write!(f, "missing prefix"),
            Self::InvalidKind(_) => write!(f, "invalid kind"),
            Self::MissingDelimiter => write!(f, "missing delimiter"),
            Self::InvalidValue(_) => write!(f, "invalid value"),
        }
    }
}

fn parse_record(mut src: &[u8]) -> Result<Record, ParseError> {
    consume_prefix(&mut src)?;
    let kind = parse_kind(&mut src).map_err(ParseError::InvalidKind)?;
    consume_delimiter(&mut src)?;
    parse_value(&mut src, kind).map_err(ParseError::InvalidValue)
}

fn consume_prefix(src: &mut &[u8]) -> Result<(), ParseError> {
    const PREFIX: u8 = b'@';

    if let Some((b, rest)) = src.split_first() {
        if *b == PREFIX {
            *src = rest;
            return Ok(());
        }
    }

    Err(ParseError::MissingPrefix)
}

fn consume_delimiter(src: &mut &[u8]) -> Result<(), ParseError> {
    const DELIMITER: u8 = b'\t';

    if let Some((b, rest)) = src.split_first() {
        if *b == DELIMITER {
            *src = rest;
            return Ok(());
        }
    }

    Err(ParseError::MissingDelimiter)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_consume_prefix() {
        let mut src = &b"@HD"[..];
        assert!(consume_prefix(&mut src).is_ok());
        assert_eq!(src, b"HD");

        let mut src = &b""[..];
        assert_eq!(consume_prefix(&mut src), Err(ParseError::MissingPrefix));

        let mut src = &b"#"[..];
        assert_eq!(consume_prefix(&mut src), Err(ParseError::MissingPrefix));
    }
}
