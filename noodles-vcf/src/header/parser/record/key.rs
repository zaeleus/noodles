use std::{error, fmt, str};

use crate::header::record::Key;

/// An error returned when a VCF header record key fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input contains invalid UTF-8.
    InvalidUtf8(str::Utf8Error),
    /// The delimiter is missing.
    MissingDelimiter,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            ParseError::InvalidUtf8(e) => Some(e),
            ParseError::MissingDelimiter => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidUtf8(_) => write!(f, "invalid UTF-8"),
            Self::MissingDelimiter => write!(f, "missing delimiter"),
        }
    }
}

pub(super) fn parse_key(src: &mut &[u8]) -> Result<Key, ParseError> {
    use memchr::memchr;

    const DELIMITER: u8 = b'=';

    let Some(i) = memchr(DELIMITER, src) else {
        return Err(ParseError::MissingDelimiter);
    };

    let (raw_key, rest) = src.split_at(i);

    let key = str::from_utf8(raw_key)
        .map(Key::from)
        .map_err(ParseError::InvalidUtf8)?;

    *src = &rest[1..];

    Ok(key)
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
        assert!(matches!(
            parse_key(&mut src),
            Err(ParseError::InvalidUtf8(_))
        ));

        let mut src = &b"fileformat"[..];
        assert_eq!(parse_key(&mut src), Err(ParseError::MissingDelimiter));
    }
}
