use std::{error, fmt};

use crate::record::QualityScore;

/// An error when a raw VCF record quality score fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
    /// The value is negative.
    Negative,
}
impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::Invalid => write!(f, "invalid input"),
            Self::Negative => write!(f, "negative value"),
        }
    }
}

pub(super) fn parse_quality_score(s: &str) -> Result<QualityScore, ParseError> {
    use crate::record::quality_score::TryFromFloatError;

    if s.is_empty() {
        return Err(ParseError::Empty);
    }

    let n = s.parse::<f32>().map_err(|_| ParseError::Invalid)?;

    QualityScore::try_from(n).map_err(|e| match e {
        TryFromFloatError::Negative => ParseError::Negative,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_position() -> Result<(), crate::record::quality_score::TryFromFloatError> {
        assert_eq!(parse_quality_score("0"), Ok(QualityScore::try_from(0.0)?));
        assert_eq!(parse_quality_score("1.0"), Ok(QualityScore::try_from(1.0)?));

        assert_eq!(parse_quality_score(""), Err(ParseError::Empty));
        assert_eq!(parse_quality_score("."), Err(ParseError::Invalid));
        assert_eq!(parse_quality_score("ndls"), Err(ParseError::Invalid));
        assert_eq!(parse_quality_score("-1.0"), Err(ParseError::Negative));

        Ok(())
    }
}
