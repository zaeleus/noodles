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
    pub fn get(&self) -> &sam::alignment::record::QualityScores {
        self.cell.get_or_init(|| get_quality_scores(&self.buf))
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

impl From<QualityScores> for sam::alignment::record::QualityScores {
    fn from(bam_quality_scores: QualityScores) -> Self {
        if let Some(sam_quality_scores) = bam_quality_scores.cell.into_inner() {
            sam_quality_scores
        } else {
            let buf = bam_quality_scores.buf;
            get_quality_scores(&buf)
        }
    }
}

fn get_quality_scores(buf: &[u8]) -> sam::alignment::record::QualityScores {
    use sam::alignment::record::quality_scores::Score;

    let scores: Vec<_> = buf
        .iter()
        .copied()
        .map(|n| Score::try_from(n).unwrap())
        .collect();

    sam::alignment::record::QualityScores::from(scores)
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
    fn test_from_quality_scores_for_sam_alignment_record_quality_scores(
    ) -> Result<(), sam::alignment::record::quality_scores::ParseError> {
        let sam_quality_scores: sam::alignment::record::QualityScores = "NDLS".parse()?;

        let bam_quality_scores = QualityScores {
            buf: BytesMut::new(),
            cell: OnceCell::from(sam_quality_scores.clone()),
        };
        let actual = sam::alignment::record::QualityScores::from(bam_quality_scores);
        assert_eq!(actual, sam_quality_scores);

        let bam_quality_scores = QualityScores {
            buf: BytesMut::from(&[45, 35, 43, 50][..]),
            cell: OnceCell::new(),
        };
        let actual = sam::alignment::record::QualityScores::from(bam_quality_scores);
        assert_eq!(actual, sam_quality_scores);

        Ok(())
    }
}
