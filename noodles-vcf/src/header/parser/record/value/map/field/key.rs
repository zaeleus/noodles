use std::{error, fmt, str};

/// An error returned when a VCF header record map value key fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input contains invalid UTF-8.
    InvalidUtf8(str::Utf8Error),
    /// Unexpected EOF
    UnexpectedEof,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidUtf8(e) => Some(e),
            Self::UnexpectedEof => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidUtf8(_) => write!(f, "invalid UTF-8"),
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
        }
    }
}

pub fn parse_key<'a>(src: &mut &'a [u8]) -> Result<&'a str, ParseError> {
    const DELIMITER: u8 = b'=';

    if let Some(i) = src.iter().position(|&b| b == DELIMITER) {
        let (buf, rest) = src.split_at(i);
        let key = str::from_utf8(buf).map_err(ParseError::InvalidUtf8)?;
        *src = &rest[1..];
        Ok(key)
    } else {
        Err(ParseError::UnexpectedEof)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_key() {
        let mut src = &b"ID="[..];
        assert_eq!(parse_key(&mut src), Ok("ID"));
        assert!(src.is_empty());

        let mut src = &b"ID"[..];
        assert_eq!(parse_key(&mut src), Err(ParseError::UnexpectedEof));
    }
}
