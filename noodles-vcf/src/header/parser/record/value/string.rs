pub mod file_format;

use std::{error, fmt, str};

pub use self::file_format::parse_file_format;

/// An error returned when a VCF header record string value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input contains invalid UTF-8.
    InvalidUtf8(str::Utf8Error),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidUtf8(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidUtf8(_) => write!(f, "invalid UTF-8"),
        }
    }
}

pub(super) fn parse_string(src: &mut &[u8]) -> Result<String, ParseError> {
    let (s, rest) = src.split_at(src.len());

    *src = rest;

    str::from_utf8(s)
        .map(String::from)
        .map_err(ParseError::InvalidUtf8)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_string() {
        let mut src = &b"noodles"[..];
        assert_eq!(parse_string(&mut src), Ok(String::from("noodles")));

        let mut src = &[0x00, 0x9f, 0x8d, 0x9c][..];
        assert!(matches!(
            parse_string(&mut src),
            Err(ParseError::InvalidUtf8(_))
        ));
    }
}
