pub mod key;
pub mod value;

use std::{borrow::Cow, error, fmt};

pub use self::{key::parse_key, value::parse_value};

/// An error returned when a VCF header record field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    UnexpectedEof,
    InvalidKey(key::ParseError),
    InvalidValue(String, value::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::UnexpectedEof => None,
            Self::InvalidKey(e) => Some(e),
            Self::InvalidValue(_, e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidKey(_) => write!(f, "invalid key"),
            Self::InvalidValue(tag, _) => write!(f, "invalid value for {tag}"),
        }
    }
}

pub fn split_field<'a>(src: &mut &'a [u8]) -> Result<Option<(&'a str, Cow<'a, str>)>, ParseError> {
    const TERMINATOR: u8 = b'>';

    if src.first().map(|&b| b == TERMINATOR).unwrap_or_default() {
        return Ok(None);
    }

    let raw_key = parse_key(src).map_err(ParseError::InvalidKey)?;
    let raw_value = parse_value(src).map_err(|e| ParseError::InvalidValue(raw_key.into(), e))?;
    consume_separator(src)?;

    Ok(Some((raw_key, raw_value)))
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
