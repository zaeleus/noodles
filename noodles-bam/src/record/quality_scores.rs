//! BAM record quality scores and iterator.
mod chars;
mod scores;

pub use self::{chars::Chars, scores::Scores};

use std::slice;

use noodles_sam::{self as sam, record::quality_scores::Score};

/// BAM record quality scores.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct QualityScores(Vec<u8>);

impl QualityScores {
    /// Creates quality scores from raw quality scores data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    /// let quality_scores = QualityScores::from(vec![45, 35, 43, 50]); // NDLS
    /// ```
    #[deprecated(
        since = "0.8.0",
        note = "Use `QualityScores::from::<Vec<u8>>` instead."
    )]
    pub fn new(qual: Vec<u8>) -> Self {
        Self::from(qual)
    }

    /// Returns the number of quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    /// let quality_scores = QualityScores::default();
    /// assert_eq!(quality_scores.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns whether there are any quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    /// let quality_scores = QualityScores::default();
    /// assert!(quality_scores.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns an iterator over quality scores as offset printable ASCII characters.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    ///
    /// let quality_scores = QualityScores::from(vec![45, 35, 43, 50]); // NDLS
    ///
    /// let mut chars = quality_scores.chars();
    ///
    /// assert_eq!(chars.next(), Some('N'));
    /// assert_eq!(chars.next(), Some('D'));
    /// assert_eq!(chars.next(), Some('L'));
    /// assert_eq!(chars.next(), Some('S'));
    /// assert_eq!(chars.next(), None);
    /// ```
    pub fn chars(&self) -> Chars<slice::Iter<'_, u8>> {
        Chars::new(self.0.iter())
    }

    /// Returns the quality score as a [`Score`] at the given index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    /// use noodles_sam::record::quality_scores::Score;
    /// let quality_scores = QualityScores::from(vec![45, 35, 43, 50]); // NDLS
    /// assert_eq!(quality_scores.get(1), Some(Score::try_from(35)));
    /// ```
    pub fn get(
        &self,
        i: usize,
    ) -> Option<Result<Score, sam::record::quality_scores::score::TryFromUByteError>> {
        self.0.get(i).copied().map(Score::try_from)
    }

    /// Returns an iterator over quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    /// use noodles_sam::record::quality_scores::Score;
    ///
    /// let quality_scores = QualityScores::from(vec![45, 35, 43, 50]); // NDLS
    ///
    /// let mut scores = quality_scores.scores();
    ///
    /// assert_eq!(scores.next(), Some(Score::try_from(45))); // N
    /// assert_eq!(scores.next(), Some(Score::try_from(35))); // D
    /// assert_eq!(scores.next(), Some(Score::try_from(43))); // L
    /// assert_eq!(scores.next(), Some(Score::try_from(50))); // S
    /// assert_eq!(scores.next(), None);
    /// ```
    pub fn scores(&self) -> Scores<slice::Iter<'_, u8>> {
        Scores::new(self.0.iter())
    }

    /// Appends a score to the end of the quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    /// use noodles_sam::record::quality_scores::Score;
    ///
    /// let mut quality_scores = QualityScores::default();
    /// quality_scores.push(Score::try_from(8)?);
    ///
    /// assert_eq!(quality_scores.as_ref(), [8]);
    /// # Ok::<_, noodles_sam::record::quality_scores::score::TryFromUByteError>(())
    /// ```
    pub fn push(&mut self, score: Score) {
        self.0.push(u8::from(score));
    }
}

impl AsRef<[u8]> for QualityScores {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<Vec<u8>> for QualityScores {
    fn as_mut(&mut self) -> &mut Vec<u8> {
        &mut self.0
    }
}

impl From<Vec<u8>> for QualityScores {
    fn from(qual: Vec<u8>) -> Self {
        Self(qual)
    }
}

impl From<Vec<Score>> for QualityScores {
    fn from(scores: Vec<Score>) -> Self {
        Self::from(scores.into_iter().map(u8::from).collect::<Vec<_>>())
    }
}

impl TryFrom<&QualityScores> for sam::record::QualityScores {
    type Error = sam::record::quality_scores::score::TryFromUByteError;

    fn try_from(quality_scores: &QualityScores) -> Result<Self, Self::Error> {
        let mut scores = Vec::with_capacity(quality_scores.len());

        for result in quality_scores.scores() {
            let score = result?;
            scores.push(score);
        }

        Ok(Self::from(scores))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_vec_score_for_quality_scores(
    ) -> Result<(), sam::record::quality_scores::score::TryFromCharError> {
        let actual = QualityScores::from(vec![
            Score::try_from('N')?,
            Score::try_from('D')?,
            Score::try_from('L')?,
            Score::try_from('S')?,
        ]);

        let expected = QualityScores::from(vec![45, 35, 43, 50]);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_try_from_quality_scores_for_sam_record_quality_scores(
    ) -> Result<(), sam::record::quality_scores::score::TryFromUByteError> {
        let quality_scores = QualityScores::from(vec![45, 35, 43, 50]); // NDLS

        let actual = sam::record::QualityScores::try_from(&quality_scores)?;
        let expected = sam::record::QualityScores::from(vec![
            Score::try_from(45)?,
            Score::try_from(35)?,
            Score::try_from(43)?,
            Score::try_from(50)?,
        ]);

        assert_eq!(actual, expected);

        Ok(())
    }
}
