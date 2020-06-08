mod score;

pub use self::score::Score;

use std::{convert::TryFrom, error, fmt, str::FromStr};

use super::NULL_FIELD;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct QualityScores {
    scores: Vec<Score>,
}

impl QualityScores {
    fn new(scores: Vec<Score>) -> Self {
        Self { scores }
    }

    pub fn scores(&self) -> &[Score] {
        &self.scores
    }

    pub fn is_empty(&self) -> bool {
        self.scores.is_empty()
    }

    pub fn len(&self) -> usize {
        self.scores.len()
    }
}

impl fmt::Display for QualityScores {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.scores.is_empty() {
            write!(f, "{}", NULL_FIELD)
        } else {
            for score in &self.scores {
                write!(f, "{}", score)?;
            }

            Ok(())
        }
    }
}

#[derive(Debug)]
pub struct ParseError(score::TryFromCharError);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid sequence: {}", self.0)
    }
}

impl FromStr for QualityScores {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.chars()
            .map(Score::try_from)
            .collect::<Result<_, _>>()
            .map(QualityScores::new)
            .map_err(ParseError)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        let scores = QualityScores::default();
        assert!(scores.is_empty());
    }

    #[test]
    fn test_len() -> Result<(), score::TryFromUByteError> {
        let sequence = [45, 35, 43, 50, 0]
            .iter()
            .cloned()
            .map(Score::try_from)
            .collect::<Result<_, _>>()
            .map(QualityScores::new)?;

        assert_eq!(sequence.len(), 5);

        Ok(())
    }

    #[test]
    fn test_fmt() -> Result<(), score::TryFromUByteError> {
        let quality_scores = [45, 35, 43, 50, 0]
            .iter()
            .cloned()
            .map(Score::try_from)
            .collect::<Result<_, _>>()
            .map(QualityScores::new)?;

        assert_eq!(quality_scores.to_string(), "NDLS!");

        Ok(())
    }

    #[test]
    fn test_from_str() -> Result<(), Box<dyn error::Error>> {
        let quality_scores = "NDLS!".parse::<QualityScores>()?;

        let expected: Vec<_> = [45, 35, 43, 50, 0]
            .iter()
            .cloned()
            .map(Score::try_from)
            .collect::<Result<_, _>>()?;

        assert_eq!(quality_scores.scores(), &expected[..]);

        Ok(())
    }
}
