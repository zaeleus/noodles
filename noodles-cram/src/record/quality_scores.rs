use std::ops::{Index, IndexMut};

use noodles_core::position::SequenceIndex;

/// Raw CRAM record quality scores.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct QualityScores(Vec<u8>);

impl QualityScores {
    /// Returns whether there are any scores.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
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
    fn from(values: Vec<u8>) -> Self {
        Self(values)
    }
}

impl<I> Index<I> for QualityScores
where
    I: SequenceIndex<u8>,
{
    type Output = I::Output;

    fn index(&self, index: I) -> &Self::Output {
        index.index(&self.0)
    }
}

impl<I> IndexMut<I> for QualityScores
where
    I: SequenceIndex<u8>,
{
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        index.index_mut(&mut self.0)
    }
}

impl From<QualityScores> for Vec<u8> {
    fn from(quality_scores: QualityScores) -> Self {
        quality_scores.0
    }
}
