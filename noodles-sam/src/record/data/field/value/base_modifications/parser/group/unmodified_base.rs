use std::{error, fmt};

use crate::record::data::field::value::base_modifications::group::UnmodifiedBase;

/// An error returned when a base modifications group unmodified base fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

pub(super) fn parse_unmodified_base(src: &mut &[u8]) -> Result<UnmodifiedBase, ParseError> {
    if let Some((raw_base, rest)) = src.split_first() {
        let base = UnmodifiedBase::try_from(*raw_base).map_err(|_| ParseError::Invalid)?;
        *src = rest;
        Ok(base)
    } else {
        Err(ParseError::UnexpectedEof)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_unmodified_base() {
        let mut src = &b"A"[..];
        assert_eq!(parse_unmodified_base(&mut src), Ok(UnmodifiedBase::A));

        let mut src = &b""[..];
        assert_eq!(
            parse_unmodified_base(&mut src),
            Err(ParseError::UnexpectedEof)
        );

        let mut src = &b"n"[..];
        assert_eq!(parse_unmodified_base(&mut src), Err(ParseError::Invalid));
    }
}
