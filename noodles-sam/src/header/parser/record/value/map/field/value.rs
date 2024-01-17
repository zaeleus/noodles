use std::{error, fmt};

use bstr::{BStr, ByteSlice};

/// An error returned when a SAM header record field tag fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is missing.
    Missing,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing => write!(f, "missing input"),
        }
    }
}

pub fn parse_value<'a>(src: &mut &'a [u8]) -> Result<&'a BStr, ParseError> {
    const DELIMITER: u8 = b'\t';

    let i = src.as_bstr().find_byte(DELIMITER).unwrap_or(src.len());
    let (buf, rest) = src.split_at(i);

    *src = rest;

    if buf.is_empty() {
        Err(ParseError::Missing)
    } else {
        Ok(buf.as_bstr())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_value() {
        let mut src = &b"ndls"[..];
        assert_eq!(parse_value(&mut src), Ok(b"ndls".as_bstr()));

        let mut src = &b""[..];
        assert_eq!(parse_value(&mut src), Err(ParseError::Missing));
    }
}
