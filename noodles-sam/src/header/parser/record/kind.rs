use std::{error, fmt};

use crate::header::record::Kind;

/// An error returned when a SAM header record kind fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

pub(super) fn parse_kind(src: &mut &[u8]) -> Result<Kind, ParseError> {
    const LEN: usize = 2;

    if src.is_empty() {
        return Err(ParseError::Empty);
    } else if src.len() < LEN {
        return Err(ParseError::Invalid);
    }

    let (raw_kind, rest) = src.split_at(LEN);
    *src = rest;

    match raw_kind {
        b"HD" => Ok(Kind::Header),
        b"SQ" => Ok(Kind::ReferenceSequence),
        b"RG" => Ok(Kind::ReadGroup),
        b"PG" => Ok(Kind::Program),
        b"CO" => Ok(Kind::Comment),
        _ => Err(ParseError::Invalid),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_kind() {
        fn t(mut src: &[u8], expected: Kind) {
            assert_eq!(parse_kind(&mut src), Ok(expected));
        }

        t(b"HD", Kind::Header);
        t(b"SQ", Kind::ReferenceSequence);
        t(b"RG", Kind::ReadGroup);
        t(b"PG", Kind::Program);
        t(b"CO", Kind::Comment);

        let mut src = &b""[..];
        assert_eq!(parse_kind(&mut src), Err(ParseError::Empty));

        let mut src = &b"H"[..];
        assert_eq!(parse_kind(&mut src), Err(ParseError::Invalid));

        let mut src = &b"ND"[..];
        assert_eq!(parse_kind(&mut src), Err(ParseError::Invalid));
    }
}
