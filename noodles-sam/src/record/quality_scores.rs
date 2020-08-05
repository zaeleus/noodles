//! SAM record quality scores and score.

mod score;

pub use self::score::Score;

use std::{convert::TryFrom, error, fmt, str::FromStr};

use super::NULL_FIELD;

/// SAM record quality scores.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct QualityScores {
    scores: Vec<Score>,
}

impl QualityScores {
    /// Creates quality scores from a list of scores.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_sam::record::{quality_scores::Score, QualityScores};
    /// let quality_scores = QualityScores::new(vec![Score::try_from(21)?]);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn new(scores: Vec<Score>) -> Self {
        Self { scores }
    }

    /// Returns the list of scores.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_sam::record::{quality_scores::Score, QualityScores};
    ///
    /// let quality_scores = QualityScores::new(vec![Score::try_from(21)?]);
    ///
    /// let actual = quality_scores.scores();
    /// let expected = [Score::try_from(21)?];
    /// assert_eq!(actual, expected);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn scores(&self) -> &[Score] {
        &self.scores
    }

    /// Returns whether the scores list is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_sam::record::{quality_scores::Score, QualityScores};
    ///
    /// let quality_scores = QualityScores::default();
    /// assert!(quality_scores.is_empty());
    ///
    /// let quality_scores = QualityScores::new(vec![Score::try_from(21)?]);
    /// assert!(!quality_scores.is_empty());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn is_empty(&self) -> bool {
        self.scores.is_empty()
    }

    /// Returns the number of scores in the list.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_sam::record::{quality_scores::Score, QualityScores};
    ///
    /// let quality_scores = QualityScores::default();
    /// assert_eq!(quality_scores.len(), 0);
    ///
    /// let quality_scores = QualityScores::new(vec![Score::try_from(21)?]);
    /// assert_eq!(quality_scores.len(), 1);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
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
            Self::Empty => f.write_str("quality scores cannot be empty"),
            Self::InvalidScore(e) => write!(f, "{}", e),
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
                .collect::<Result<_, _>>()
                .map(QualityScores::new)
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
    fn test_from_str() -> Result<(), score::TryFromUByteError> {
        let expected = [45, 35, 43, 50, 0]
            .iter()
            .cloned()
            .map(Score::try_from)
            .collect::<Result<_, _>>()
            .map(QualityScores::new)?;
        assert_eq!("NDLS!".parse(), Ok(expected));

        assert_eq!("*".parse::<QualityScores>(), Ok(QualityScores::default()));

        assert_eq!("".parse::<QualityScores>(), Err(ParseError::Empty));

        Ok(())
    }
}
