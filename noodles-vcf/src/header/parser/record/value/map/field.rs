pub mod key;
pub mod value;

use std::{error, fmt};

pub use self::{key::parse_key, value::parse_value};

/// An error returned when a VCF header record field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    UnexpectedEof,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
        }
    }
}

pub fn consume_separator(src: &mut &[u8]) -> Result<bool, ParseError> {
    const SEPARATOR: u8 = b',';

    if let Some((b, rest)) = src.split_first() {
        if *b == SEPARATOR {
            *src = rest;
            Ok(true)
        } else {
            Ok(false)
        }
    } else {
        Err(ParseError::UnexpectedEof)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_consume_separator() {
        let mut src = &b","[..];
        assert_eq!(consume_separator(&mut src), Ok(true));
        assert!(src.is_empty());

        let mut src = &b">"[..];
        assert_eq!(consume_separator(&mut src), Ok(false));
        assert_eq!(src, b">");

        let mut src = &b""[..];
        assert_eq!(consume_separator(&mut src), Err(ParseError::UnexpectedEof));
    }
}
