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

impl<'a> sam::alignment::record::QualityScores for QualityScores<'a> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_> {
        Box::new(self.as_ref().iter().copied())
    }
}

impl<'a> AsRef<[u8]> for QualityScores<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> From<QualityScores<'a>> for sam::alignment::record_buf::QualityScores {
    fn from(quality_scores: QualityScores<'a>) -> Self {
        Self::from(quality_scores.0.to_vec())
    }
}
