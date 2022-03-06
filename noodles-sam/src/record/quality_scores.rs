//! SAM record quality scores and score.

pub mod score;

pub use self::score::Score;

use std::{error, fmt, ops::Deref, str::FromStr};

use super::NULL_FIELD;

/// SAM record quality scores.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct QualityScores(Vec<Score>);

impl QualityScores {
    /// Removes all scores.
    ///
    /// This does not affect the capacity of the list.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{quality_scores::Score, QualityScores};
    ///
    /// let mut quality_scores = QualityScores::from(vec![Score::try_from('!')?]);
    /// assert!(!quality_scores.is_empty());
    ///
    /// quality_scores.clear();
    /// assert!(quality_scores.is_empty());
    /// # Ok::<_, noodles_sam::record::quality_scores::score::TryFromCharError>(())
    /// ```
    pub fn clear(&mut self) {
        self.0.clear();
    }
}

impl Deref for QualityScores {
    type Target = [Score];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for QualityScores {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            write!(f, "{}", NULL_FIELD)
        } else {
            for score in self.iter() {
                write!(f, "{}", score)?;
            }

            Ok(())
        }
    }
}

impl From<Vec<Score>> for QualityScores {
    fn from(scores: Vec<Score>) -> Self {
        Self(scores)
    }
}

/// An error returned when raw SAM record quality scores fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The raw quality scores has an invalid score.
    InvalidScore(score::TryFromCharError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidScore(e) => write!(f, "invalid score: {}", e),
        }
    }
}

impl FromStr for QualityScores {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            NULL_FIELD => Ok(Self::default()),
            _ => s
                .chars()
                .map(Score::try_from)
                .collect::<Result<Vec<_>, _>>()
                .map(Self::from)
                .map_err(ParseError::InvalidScore),
        }
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
            .collect::<Result<Vec<_>, _>>()
            .map(QualityScores::from)?;

        assert_eq!(sequence.len(), 5);

        Ok(())
    }

    #[test]
    fn test_fmt() -> Result<(), score::TryFromUByteError> {
        let quality_scores = [45, 35, 43, 50, 0]
            .iter()
            .cloned()
            .map(Score::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map(QualityScores::from)?;

        assert_eq!(quality_scores.to_string(), "NDLS!");

        Ok(())
    }

    #[test]
    fn test_from_str() -> Result<(), score::TryFromUByteError> {
        let expected = [45, 35, 43, 50, 0]
            .iter()
            .cloned()
            .map(Score::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map(QualityScores::from)?;
        assert_eq!("NDLS!".parse(), Ok(expected));

        assert_eq!("*".parse::<QualityScores>(), Ok(QualityScores::default()));

        assert_eq!("".parse::<QualityScores>(), Err(ParseError::Empty));

        Ok(())
    }
}
