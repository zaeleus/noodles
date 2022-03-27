//! SAM record quality scores and score.

pub mod score;

pub use self::score::Score;

use std::{
    error, fmt,
    ops::{Index, IndexMut},
    str::FromStr,
};

use noodles_core::position::SequenceIndex;

/// SAM record quality scores.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct QualityScores(Vec<Score>);

impl QualityScores {
    /// Returns whether there are any scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::QualityScores;
    /// let quality_scores = QualityScores::default();
    /// assert!(quality_scores.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::QualityScores;
    /// let quality_scores = QualityScores::default();
    /// assert_eq!(quality_scores.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }

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

    /// Returns a reference to the score at the given index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam::record::{quality_scores::Score, QualityScores};
    ///
    /// let quality_scores: QualityScores = "NDLS".parse()?;
    ///
    /// let i = Position::try_from(2)?;
    /// assert_eq!(quality_scores.get(i), Some(&Score::try_from('D')?));
    ///
    /// let i = Position::try_from(8)?;
    /// assert!(quality_scores.get(i).is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: SequenceIndex<Score>,
    {
        index.get(self.0.as_ref())
    }

    /// Returns a mutable reference to the score at the given index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam::record::{quality_scores::Score, QualityScores};
    ///
    /// let mut quality_scores: QualityScores = "NDLS".parse()?;
    ///
    /// let i = Position::try_from(2)?;
    /// if let Some(score) = quality_scores.get_mut(i) {
    ///     *score = Score::try_from('!')?;
    /// }
    ///
    /// assert_eq!(quality_scores.get(i), Some(&Score::try_from('!')?));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn get_mut<I>(&mut self, index: I) -> Option<&mut I::Output>
    where
        I: SequenceIndex<Score>,
    {
        index.get_mut(self.0.as_mut())
    }

    /// Appends a score.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{quality_scores::Score, QualityScores};
    ///
    /// let mut quality_scores = QualityScores::from(vec![Score::try_from('N')?]);
    /// quality_scores.push(Score::try_from('D')?);
    ///
    /// let expected = "ND".parse()?;
    /// assert_eq!(quality_scores, expected);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn push(&mut self, score: Score) {
        self.0.push(score);
    }
}

impl AsRef<[Score]> for QualityScores {
    fn as_ref(&self) -> &[Score] {
        &self.0
    }
}

impl AsMut<Vec<Score>> for QualityScores {
    fn as_mut(&mut self) -> &mut Vec<Score> {
        &mut self.0
    }
}

impl fmt::Display for QualityScores {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for score in self.as_ref() {
            write!(f, "{}", score)?;
        }

        Ok(())
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
        if s.is_empty() {
            Err(ParseError::Empty)
        } else {
            s.chars()
                .map(Score::try_from)
                .collect::<Result<Vec<_>, _>>()
                .map(Self::from)
                .map_err(ParseError::InvalidScore)
        }
    }
}

impl<I> Index<I> for QualityScores
where
    I: SequenceIndex<Score>,
{
    type Output = I::Output;

    fn index(&self, index: I) -> &Self::Output {
        index.index(&self.0)
    }
}

impl<I> IndexMut<I> for QualityScores
where
    I: SequenceIndex<Score>,
{
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        index.index_mut(&mut self.0)
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
        let quality_scores = QualityScores::default();
        assert!(quality_scores.to_string().is_empty());

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

        assert_eq!(
            "*".parse::<QualityScores>(),
            Ok(QualityScores::from(vec![Score::try_from(9)?]))
        );

        assert_eq!("".parse::<QualityScores>(), Err(ParseError::Empty));

        Ok(())
    }
}
