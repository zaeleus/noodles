use std::{error, fmt, str};

use crate::header::record::Key;

/// An error returned when a VCF header record key fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid(str::Utf8Error),
    /// The delimiter is missing.
    MissingDelimiter,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            ParseError::Invalid(e) => Some(e),
            ParseError::MissingDelimiter => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(_) => write!(f, "invalid input"),
            Self::MissingDelimiter => write!(f, "missing delimiter"),
        }
    }
}

pub(super) fn parse_key(src: &mut &[u8]) -> Result<Key, ParseError> {
    const DELIMITER: u8 = b'=';

    if let Some(i) = src.iter().position(|&b| b == DELIMITER) {
        let (raw_key, rest) = src.split_at(i);

        let key = str::from_utf8(raw_key)
            .map(Key::from)
            .map_err(ParseError::Invalid)?;

        *src = &rest[1..];

        Ok(key)
    } else {
        Err(ParseError::MissingDelimiter)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_key() {
        use crate::header::record::key;

        let mut src = &b"fileformat="[..];
        assert_eq!(parse_key(&mut src), Ok(key::FILE_FORMAT));

        let mut src = &[0x00, 0x9f, 0x8d, 0x9c, b'='][..];
        assert!(matches!(parse_key(&mut src), Err(ParseError::Invalid(_))));

        let mut src = &b"fileformat"[..];
        assert_eq!(parse_key(&mut src), Err(ParseError::MissingDelimiter));
    }
}
