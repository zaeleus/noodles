use std::io;

use noodles_sam as sam;

/// Raw BAM record quality scores.
#[derive(Debug, Eq, PartialEq)]
pub struct QualityScores<'a>(&'a [u8]);

impl<'a> QualityScores<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }

    /// Returns whether there are any scores.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of scores.
    pub fn len(&self) -> usize {
        self.0.len()
    }
}

impl<'a> AsRef<[u8]> for QualityScores<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryFrom<QualityScores<'a>> for sam::record::QualityScores {
    type Error = io::Error;

    fn try_from(bam_quality_scores: QualityScores<'a>) -> Result<Self, Self::Error> {
        use crate::reader::record::get_quality_scores;

        let mut src = bam_quality_scores.0;
        let mut quality_scores = Self::default();
        let base_count = src.len();
        get_quality_scores(&mut src, &mut quality_scores, base_count)?;

        Ok(quality_scores)
    }
}
