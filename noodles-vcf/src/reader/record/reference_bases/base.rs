use std::{error, fmt};

use noodles_core as core;

use crate::record::reference_bases::Base;

/// An error when a raw VCF record base fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid(char),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(c) => write!(f, "expected {{A, C, G, T, N}}, got {}", c),
        }
    }
}

impl From<ParseError> for core::Error {
    fn from(e: ParseError) -> Self {
        Self::new(core::error::Kind::Parse, e)
    }
}

pub(super) fn parse_base(c: char) -> Result<Base, ParseError> {
    match c.to_ascii_uppercase() {
        'A' => Ok(Base::A),
        'C' => Ok(Base::C),
        'G' => Ok(Base::G),
        'T' => Ok(Base::T),
        'N' => Ok(Base::N),
        _ => Err(ParseError::Invalid(c)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_base() {
        assert_eq!(parse_base('A'), Ok(Base::A));
        assert_eq!(parse_base('C'), Ok(Base::C));
        assert_eq!(parse_base('G'), Ok(Base::G));
        assert_eq!(parse_base('T'), Ok(Base::T));
        assert_eq!(parse_base('N'), Ok(Base::N));

        assert_eq!(parse_base('a'), Ok(Base::A));
        assert_eq!(parse_base('c'), Ok(Base::C));
        assert_eq!(parse_base('g'), Ok(Base::G));
        assert_eq!(parse_base('t'), Ok(Base::T));
        assert_eq!(parse_base('n'), Ok(Base::N));

        assert_eq!(parse_base('Z'), Err(ParseError::Invalid('Z')));
        assert_eq!(parse_base('z'), Err(ParseError::Invalid('z')));
    }
}
