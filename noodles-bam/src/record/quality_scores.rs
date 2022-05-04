//! BAM record quality scores.

use std::str::FromStr;

use bytes::BytesMut;
use noodles_sam::{self as sam, alignment::record::AlignmentQualityScores};
use once_cell::sync::OnceCell;

/// Lazy BAM record quality scores.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct QualityScores {
    pub(crate) buf: BytesMut,
    cell: OnceCell<sam::alignment::record::QualityScores>,
}

impl QualityScores {
    /// Returns the number of scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    ///
    /// let quality_scores = QualityScores::default();
    /// assert_eq!(quality_scores.len(), 0);
    ///
    /// let quality_scores: QualityScores = "NDLS".parse()?;
    /// assert_eq!(quality_scores.len(), 4);
    /// # Ok::<_, noodles_bam::record::quality_scores::ParseError>(())
    /// ```
    pub fn len(&self) -> usize {
        self.cell
            .get()
            .map(|quality_scores| quality_scores.len())
            .unwrap_or(self.buf.len())
    }

    /// Returns whether there are any scores.
    pub fn is_empty(&self) -> bool {
        self.cell
            .get()
            .map(|quality_scores| quality_scores.is_empty())
            .unwrap_or(self.buf.is_empty())
    }

    /// Removes all scores.
    pub fn clear(&mut self) {
        self.buf = BytesMut::new();
        self.cell.take();
    }

    /// Returns alignment record quality scores.
    pub fn try_get(&self) -> Result<&sam::alignment::record::QualityScores, ParseError> {
        self.cell
            .get_or_try_init(|| sam::alignment::record::QualityScores::try_from(self.buf.to_vec()))
    }
}

/// An error returned when raw alignment record quality scores fail to parse.
pub type ParseError = sam::alignment::record::quality_scores::ParseError;

impl FromStr for QualityScores {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let sam_quality_scores: sam::alignment::record::QualityScores = s.parse()?;
        Ok(Self::from(sam_quality_scores))
    }
}

impl From<sam::alignment::record::QualityScores> for QualityScores {
    fn from(quality_scores: sam::alignment::record::QualityScores) -> Self {
        let cell = OnceCell::new();
        cell.set(quality_scores).ok();

        QualityScores {
            buf: BytesMut::new(),
            cell,
        }
    }
}

impl TryFrom<QualityScores> for sam::alignment::record::QualityScores {
    type Error = ParseError;

    fn try_from(bam_quality_scores: QualityScores) -> Result<Self, Self::Error> {
        if let Some(sam_quality_scores) = bam_quality_scores.cell.into_inner() {
            Ok(sam_quality_scores)
        } else {
            let buf = bam_quality_scores.buf;
            Self::try_from(buf.to_vec())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_sam_alignment_record_quality_scores_for_quality_scores(
    ) -> Result<(), sam::alignment::record::quality_scores::ParseError> {
        let sam_quality_scores: sam::alignment::record::QualityScores = "NDLS".parse()?;

        let actual = QualityScores::from(sam_quality_scores.clone());
        assert!(actual.buf.is_empty());
        assert_eq!(actual.cell.get(), Some(&sam_quality_scores));

        Ok(())
    }

    #[test]
    fn test_from_try_quality_scores_for_sam_alignment_record_quality_scores(
    ) -> Result<(), sam::alignment::record::quality_scores::ParseError> {
        let sam_quality_scores: sam::alignment::record::QualityScores = "NDLS".parse()?;

        let bam_quality_scores = QualityScores {
            buf: BytesMut::new(),
            cell: OnceCell::from(sam_quality_scores.clone()),
        };
        let actual = sam::alignment::record::QualityScores::try_from(bam_quality_scores)?;
        assert_eq!(actual, sam_quality_scores);

        let bam_quality_scores = QualityScores {
            buf: BytesMut::from(&[45, 35, 43, 50][..]),
            cell: OnceCell::new(),
        };
        let actual = sam::alignment::record::QualityScores::try_from(bam_quality_scores)?;
        assert_eq!(actual, sam_quality_scores);

        Ok(())
    }
}
