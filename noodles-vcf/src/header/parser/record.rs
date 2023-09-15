#![allow(dead_code)]

mod key;
mod value;

use std::{error, fmt};

use self::{key::parse_key, value::parse_value};
use crate::header::{record::Key, FileFormat, Record};

/// An error returned when a VCF header record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The header prefix is missing.
    MissingPrefix,
    /// The key is invalid.
    InvalidKey(key::ParseError),
    /// The value is invalid.
    InvalidValue(Key, value::ParseError),
}

impl ParseError {
    pub(super) fn key(&self) -> Option<&Key> {
        match self {
            Self::InvalidValue(k, _) => Some(k),
            _ => None,
        }
    }
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::MissingPrefix => None,
            Self::InvalidKey(e) => Some(e),
            Self::InvalidValue(_, e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingPrefix => write!(f, "missing prefix"),
            Self::InvalidKey(_) => write!(f, "invalid key"),
            Self::InvalidValue(..) => write!(f, "invalid value"),
        }
    }
}

#[allow(missing_docs)]
pub fn parse_record(mut src: &[u8], file_format: FileFormat) -> Result<Record, ParseError> {
    consume_prefix(&mut src)?;
    let key = parse_key(&mut src).map_err(ParseError::InvalidKey)?;
    parse_value(&mut src, file_format, key.clone()).map_err(|e| ParseError::InvalidValue(key, e))
}

fn consume_prefix(src: &mut &[u8]) -> Result<(), ParseError> {
    const PREFIX: &[u8] = b"##";

    if let Some(rest) = src.strip_prefix(PREFIX) {
        *src = rest;
        Ok(())
    } else {
        Err(ParseError::MissingPrefix)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_consume_prefix() {
        let mut src = &b"##fileformat"[..];
        assert!(consume_prefix(&mut src).is_ok());
        assert_eq!(src, b"fileformat");

        let mut src = &b""[..];
        assert_eq!(consume_prefix(&mut src), Err(ParseError::MissingPrefix));

        let mut src = &b"#"[..];
        assert_eq!(consume_prefix(&mut src), Err(ParseError::MissingPrefix));

        let mut src = &b"@"[..];
        assert_eq!(consume_prefix(&mut src), Err(ParseError::MissingPrefix));
    }
}
