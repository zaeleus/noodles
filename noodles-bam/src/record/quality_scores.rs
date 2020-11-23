//! BAM record quality scores and iterator.
mod chars;

pub use self::chars::Chars;

use std::{ops::Deref, slice};

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
    ///
    /// // NDLS
    /// let data = [45, 35, 43, 50];
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
    /// // NDLS
    /// let data = [45, 35, 43, 50];
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
}

impl<'a> Deref for QualityScores<'a> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.qual
    }
}
