//! VCF record quality score.

use std::{convert::TryFrom, error, fmt, num, ops::Deref, str::FromStr};

use super::{value::parse_f32_case_insensitive_extended, MISSING_FIELD};

const MIN: f32 = 0.0;

/// A VCF record quality score.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct QualityScore(Option<f32>);

impl Deref for QualityScore {
    type Target = Option<f32>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for QualityScore {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            Some(score) => write!(f, "{}", score),
            None => f.write_str(MISSING_FIELD),
        }
    }
}

/// An error returned when a raw float fails to convert to a quality score.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromFloatError {
    /// The input is negative.
    Negative,
}

impl error::Error for TryFromFloatError {}

impl fmt::Display for TryFromFloatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Negative => f.write_str("negative value"),
        }
    }
}

impl TryFrom<f32> for QualityScore {
    type Error = TryFromFloatError;

    fn try_from(value: f32) -> Result<Self, Self::Error> {
        if value < MIN {
            Err(TryFromFloatError::Negative)
        } else {
            Ok(Self(Some(value)))
        }
    }
}

/// An error returned when a raw VCF record quality score cannot be parsed.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid(num::ParseFloatError),
    /// The quality score is invalid.
    InvalidValue(TryFromFloatError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid(e) => write!(f, "invalid input: {}", e),
            Self::InvalidValue(e) => write!(f, "invalid value: {}", e),
        }
    }
}

impl FromStr for QualityScore {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            MISSING_FIELD => Ok(Self(None)),
            _ => parse_f32_case_insensitive_extended(s)
                .map_err(ParseError::Invalid)
                .and_then(|value| Self::try_from(value).map_err(ParseError::InvalidValue)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let quality_score = QualityScore(Some(1.3));
        assert_eq!(quality_score.to_string(), "1.3");

        let quality_score = QualityScore(None);
        assert_eq!(quality_score.to_string(), ".");
    }

    #[test]
    fn test_try_from_float_for_quality_score() {
        assert_eq!(QualityScore::try_from(0.0), Ok(QualityScore(Some(0.0))));
        assert_eq!(QualityScore::try_from(13.0), Ok(QualityScore(Some(13.0))));

        assert_eq!(
            QualityScore::try_from(-8.0),
            Err(TryFromFloatError::Negative)
        );
    }

    #[test]
    fn test_from_str() {
        assert_eq!(".".parse(), Ok(QualityScore(None)));
        assert_eq!("5.8".parse(), Ok(QualityScore(Some(5.8))));
        assert_eq!("Infinity".parse(), Ok(QualityScore(Some(f32::INFINITY))));

        assert_eq!("".parse::<QualityScore>(), Err(ParseError::Empty));
        assert!(matches!(
            "ndls".parse::<QualityScore>(),
            Err(ParseError::Invalid(_))
        ));
        assert!(matches!(
            "-8.5".parse::<QualityScore>(),
            Err(ParseError::InvalidValue(_))
        ));
    }
}
