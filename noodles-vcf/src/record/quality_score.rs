//! VCF record quality score.

use std::{error, fmt, num, str::FromStr};

const MIN: f32 = 0.0;

/// A VCF record quality score.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct QualityScore(f32);

impl fmt::Display for QualityScore {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
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
            Ok(Self(value))
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

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Empty => None,
            Self::Invalid(e) => Some(e),
            Self::InvalidValue(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid(_) => f.write_str("invalid input"),
            Self::InvalidValue(_) => f.write_str("invalid value"),
        }
    }
}

impl FromStr for QualityScore {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            _ => s
                .parse::<f32>()
                .map_err(ParseError::Invalid)
                .and_then(|value| Self::try_from(value).map_err(ParseError::InvalidValue)),
        }
    }
}

impl From<QualityScore> for f32 {
    fn from(quality_score: QualityScore) -> Self {
        quality_score.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let quality_score = QualityScore(1.3);
        assert_eq!(quality_score.to_string(), "1.3");
    }

    #[test]
    fn test_try_from_float_for_quality_score() {
        assert_eq!(QualityScore::try_from(0.0), Ok(QualityScore(0.0)));
        assert_eq!(QualityScore::try_from(13.0), Ok(QualityScore(13.0)));

        assert_eq!(
            QualityScore::try_from(-8.0),
            Err(TryFromFloatError::Negative)
        );
    }

    #[test]
    fn test_from_str() {
        assert_eq!("5.8".parse(), Ok(QualityScore(5.8)));
        assert_eq!("Infinity".parse(), Ok(QualityScore(f32::INFINITY)));

        assert_eq!("".parse::<QualityScore>(), Err(ParseError::Empty));
        assert!(matches!(
            ".".parse::<QualityScore>(),
            Err(ParseError::Invalid(_))
        ));
        assert!(matches!(
            "ndls".parse::<QualityScore>(),
            Err(ParseError::Invalid(_))
        ));
        assert!(matches!(
            "-8.5".parse::<QualityScore>(),
            Err(ParseError::InvalidValue(_))
        ));
    }

    #[test]
    fn test_from_quality_score_for_f32() {
        assert_eq!(f32::from(QualityScore(13.0)).to_bits(), (13.0f32).to_bits());
    }
}
