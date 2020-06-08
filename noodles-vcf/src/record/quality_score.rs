use std::{error, fmt, ops::Deref, str::FromStr};

use super::MISSING_FIELD;

#[derive(Clone, Debug, Default, PartialEq)]
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

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid quality score: {}", self.0)
    }
}

impl FromStr for QualityScore {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError(s.into())),
            MISSING_FIELD => Ok(Self(None)),
            _ => s
                .parse()
                .map(|value| Self(Some(value)))
                .map_err(|_| ParseError(s.into())),
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
    fn test_from_str() -> Result<(), ParseError> {
        assert!(".".parse::<QualityScore>()?.is_none());

        assert_eq!("5.8".parse::<QualityScore>()?, QualityScore(Some(5.8)));

        assert!("".parse::<QualityScore>().is_err());
        assert!("ndls".parse::<QualityScore>().is_err());

        Ok(())
    }
}
