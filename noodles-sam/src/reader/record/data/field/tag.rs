use std::{error, fmt};

use crate::record::data::field::{tag, Tag};

/// An error when a raw BAM record data field tag fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The input is invalid.
    Invalid(tag::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::UnexpectedEof => None,
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::Invalid(_) => write!(f, "invalid input"),
        }
    }
}

pub(super) fn parse_tag(src: &mut &[u8]) -> Result<Tag, ParseError> {
    if src.len() < 2 {
        return Err(ParseError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(2);
    *src = rest;

    Tag::try_from(buf).map_err(ParseError::Invalid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tag() {
        let mut src = &b"NH"[..];
        assert_eq!(parse_tag(&mut src), Ok(tag::ALIGNMENT_HIT_COUNT));

        let mut src = &b""[..];
        assert_eq!(parse_tag(&mut src), Err(ParseError::UnexpectedEof));

        let mut src = &b"X!"[..];
        assert!(matches!(parse_tag(&mut src), Err(ParseError::Invalid(_))));
    }
}
