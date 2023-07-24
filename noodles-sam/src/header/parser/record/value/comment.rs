use std::{error, fmt, str};

/// An error returned when a SAM header record comment value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The delimiter is invalid.
    InvalidDelimiter,
    /// The input is invalid.
    Invalid(str::Utf8Error),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidDelimiter => None,
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDelimiter => write!(f, "invalid delimiter"),
            Self::Invalid(_) => write!(f, "invalid input"),
        }
    }
}

pub(super) fn parse_comment(src: &mut &[u8]) -> Result<String, ParseError> {
    consume_delimiter(src)?;

    let (buf, rest) = src.split_at(src.len());

    let s = str::from_utf8(buf)
        .map(String::from)
        .map_err(ParseError::Invalid)?;

    *src = rest;

    Ok(s)
}

fn consume_delimiter(src: &mut &[u8]) -> Result<(), ParseError> {
    const PREFIX: u8 = b'\t';

    if let Some((b, rest)) = src.split_first() {
        if *b == PREFIX {
            *src = rest;
            return Ok(());
        }
    }

    Err(ParseError::InvalidDelimiter)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_comment() {
        let mut src = &b"\tnoodles"[..];
        assert_eq!(parse_comment(&mut src), Ok(String::from("noodles")));

        let mut src = &b"noodles"[..];
        assert_eq!(parse_comment(&mut src), Err(ParseError::InvalidDelimiter));

        let mut src = &[b'\t', 0x00, 0x9f, 0x8d, 0x9c][..];
        assert!(matches!(
            parse_comment(&mut src),
            Err(ParseError::Invalid(_))
        ));
    }
}
