use std::{error, fmt};

/// An error when a raw VCF record quality score fails to parse.
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

pub(super) fn parse_quality_score(s: &str) -> Result<f32, ParseError> {
    if s.is_empty() {
        return Err(ParseError::Empty);
    }

    s.parse().map_err(|_| ParseError::Invalid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_position() {
        assert_eq!(parse_quality_score("0"), Ok(0.0));
        assert_eq!(parse_quality_score("1.0"), Ok(1.0));

        assert_eq!(parse_quality_score(""), Err(ParseError::Empty));
        assert_eq!(parse_quality_score("."), Err(ParseError::Invalid));
        assert_eq!(parse_quality_score("ndls"), Err(ParseError::Invalid));
    }
}
