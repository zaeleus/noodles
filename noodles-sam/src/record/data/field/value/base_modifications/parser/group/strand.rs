use std::{error, fmt};

use crate::record::data::field::value::base_modifications::group::Strand;

/// An error returned when a base modifications group strand fails to parse.
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

pub fn parse_strand(src: &mut &[u8]) -> Result<Strand, ParseError> {
    if let Some((raw_strand, rest)) = src.split_first() {
        let strand = Strand::try_from(*raw_strand).map_err(|_| ParseError::Invalid)?;
        *src = rest;
        Ok(strand)
    } else {
        Err(ParseError::UnexpectedEof)
    }
}

#[cfg(test)]
mod tests {
    use crate::record::data::field::value::base_modifications::group::Strand;

    use super::*;

    #[test]
    fn test_parse_strand() {
        let mut src = &b"+"[..];
        assert_eq!(parse_strand(&mut src), Ok(Strand::Forward));

        let mut src = &b"-"[..];
        assert_eq!(parse_strand(&mut src), Ok(Strand::Reverse));

        let mut src = &b""[..];
        assert_eq!(parse_strand(&mut src), Err(ParseError::UnexpectedEof));

        let mut src = &b"."[..];
        assert_eq!(parse_strand(&mut src), Err(ParseError::Invalid));
    }
}
