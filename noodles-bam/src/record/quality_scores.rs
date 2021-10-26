//! BAM record quality scores and iterator.
mod chars;
mod scores;

pub use self::{chars::Chars, scores::Scores};

use std::{ops::Deref, slice};

use noodles_sam as sam;

/// BAM record quality scores.
#[derive(Debug)]
pub struct QualityScores<'a> {
    qual: &'a [u8],
}

impl<'a> QualityScores<'a> {
    /// Creates quality scores from raw quality scores data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    /// let data = [45, 35, 43, 50]; // NDLS
    /// let quality_scores = QualityScores::new(&data);
    /// ```
    pub fn new(qual: &'a [u8]) -> Self {
        Self { qual }
    }

    /// Returns an iterator over quality scores as offset printable ASCII characters.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    ///
    /// let data = [45, 35, 43, 50]; // NDLS
    /// let quality_scores = QualityScores::new(&data);
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
        Chars::new(self.qual.iter())
    }

    /// Returns an iterator over quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    /// use noodles_sam::record::quality_scores::Score;
    ///
    /// let data = [45, 35, 43, 50]; // NDLS
    /// let quality_scores = QualityScores::new(&data);
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
        Scores::new(self.qual.iter())
    }
}

impl<'a> Deref for QualityScores<'a> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.qual
    }
}

impl<'a> TryFrom<QualityScores<'a>> for sam::record::QualityScores {
    type Error = sam::record::quality_scores::score::TryFromUByteError;

    fn try_from(quality_scores: QualityScores<'_>) -> Result<Self, Self::Error> {
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

        let data = [45, 35, 43, 50]; // NDLS
        let quality_scores = QualityScores::new(&data);

        let actual = sam::record::QualityScores::try_from(quality_scores)?;
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
