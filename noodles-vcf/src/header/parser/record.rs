#![allow(dead_code)]

mod key;

use std::{error, fmt};

use self::key::parse_key;
use crate::header::Record;

/// An error returned when a VCF header record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The header prefix is missing.
    MissingPrefix,
    /// The key is invalid.
    InvalidKey(key::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::MissingPrefix => None,
            Self::InvalidKey(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingPrefix => write!(f, "missing prefix"),
            Self::InvalidKey(_) => write!(f, "invalid key"),
        }
    }
}

pub(super) fn parse_record(mut src: &[u8]) -> Result<Record, ParseError> {
    consume_prefix(&mut src)?;
    let _ = parse_key(&mut src).map_err(ParseError::InvalidKey)?;
    todo!()
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
