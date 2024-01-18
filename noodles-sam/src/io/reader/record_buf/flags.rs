use std::{error, fmt};

use crate::alignment::record::Flags;

/// An error when raw SAM record flags fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    Invalid(lexical_core::Error),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(_) => write!(f, "invalid input"),
        }
    }
}

pub(crate) fn parse_flags(src: &[u8]) -> Result<Flags, ParseError> {
    lexical_core::parse::<u16>(src)
        .map_err(ParseError::Invalid)
        .map(Flags::from)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_flags() {
        assert_eq!(parse_flags(b"0"), Ok(Flags::empty()));
        assert_eq!(parse_flags(b"4"), Ok(Flags::UNMAPPED));

        assert!(matches!(parse_flags(b""), Err(ParseError::Invalid(_))));
        assert!(matches!(parse_flags(b"-4"), Err(ParseError::Invalid(_))));
        assert!(matches!(parse_flags(b"n"), Err(ParseError::Invalid(_))));
    }
}
