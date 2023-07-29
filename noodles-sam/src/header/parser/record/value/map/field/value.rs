use std::{error, fmt, str};

/// An error returned when a SAM header record field tag fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is missing.
    Missing,
    /// The input is invalid.
    Invalid(str::Utf8Error),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Missing => None,
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing => write!(f, "missing input"),
            Self::Invalid(_) => write!(f, "invalid input"),
        }
    }
}

pub fn parse_value<'a>(src: &mut &'a [u8]) -> Result<&'a str, ParseError> {
    use memchr::memchr;

    const DELIMITER: u8 = b'\t';

    let i = memchr(DELIMITER, src).unwrap_or(src.len());
    let (buf, rest) = src.split_at(i);

    *src = rest;

    if buf.is_empty() {
        Err(ParseError::Missing)
    } else {
        str::from_utf8(buf).map_err(ParseError::Invalid)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_value() {
        let mut src = &b"ndls"[..];
        assert_eq!(parse_value(&mut src), Ok("ndls"));

        let mut src = &b""[..];
        assert_eq!(parse_value(&mut src), Err(ParseError::Missing));

        let mut src = &[0x00, 0x9f, 0x8d, 0x9c][..];
        assert!(matches!(parse_value(&mut src), Err(ParseError::Invalid(_))));
    }
}
