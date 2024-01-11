use std::{error, fmt};

use noodles_core::Position;

/// An error when a raw SAM record position fail to parse.
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

pub(super) fn parse_alignment_start(src: &[u8]) -> Result<Option<Position>, ParseError> {
    lexical_core::parse(src)
        .map_err(ParseError::Invalid)
        .map(Position::new)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_alignment_start() -> Result<(), noodles_core::position::TryFromIntError> {
        assert_eq!(parse_alignment_start(b"0"), Ok(None));
        assert_eq!(
            parse_alignment_start(b"8"),
            Ok(Some(Position::try_from(8)?))
        );

        assert!(matches!(
            parse_alignment_start(b"-1"),
            Err(ParseError::Invalid(_))
        ));
        assert!(matches!(
            parse_alignment_start(b"n"),
            Err(ParseError::Invalid(_))
        ));

        Ok(())
    }
}
