//! BAM record quality scores and iterator.
mod chars;
mod scores;

pub use self::{chars::Chars, scores::Scores};

use std::{
    ops::{Deref, DerefMut},
    slice,
};

use noodles_sam as sam;

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
    /// let data = vec![45, 35, 43, 50]; // NDLS
    /// let quality_scores = QualityScores::from(data);
    /// ```
    #[deprecated(
        since = "0.8.0",
        note = "Use `QualityScores::from::<Vec<u8>>` instead."
    )]
    pub fn new(qual: Vec<u8>) -> Self {
        Self::from(qual)
    }

    /// Returns an iterator over quality scores as offset printable ASCII characters.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    ///
    /// let data = vec![45, 35, 43, 50]; // NDLS
    /// let quality_scores = QualityScores::from(data);
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
        Chars::new(self.iter())
    }

    /// Returns an iterator over quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    /// use noodles_sam::record::quality_scores::Score;
    ///
    /// let data = vec![45, 35, 43, 50]; // NDLS
    /// let quality_scores = QualityScores::from(data);
    ///
    /// let mut scores = quality_scores.scores();
    ///
    /// assert_eq!(scores.next(), Some(Score::try_from(45)));
    /// assert_eq!(scores.next(), Some(Score::try_from(35)));
    /// assert_eq!(scores.next(), Some(Score::try_from(43)));
    /// assert_eq!(scores.next(), Some(Score::try_from(50)));
    /// assert_eq!(scores.next(), None);
    /// ```
    pub fn scores(&self) -> Scores<slice::Iter<'_, u8>> {
        Scores::new(self.iter())
    }
}

impl Deref for QualityScores {
    type Target = Vec<u8>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for QualityScores {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<Vec<u8>> for QualityScores {
    fn from(qual: Vec<u8>) -> Self {
        Self(qual)
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
    fn test_try_from_quality_scores_for_sam_record_quality_scores(
    ) -> Result<(), sam::record::quality_scores::score::TryFromUByteError> {
        use sam::record::quality_scores::Score;

        let data = vec![45, 35, 43, 50]; // NDLS
        let quality_scores = QualityScores::from(data);

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
