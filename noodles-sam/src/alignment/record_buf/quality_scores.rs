use std::{
    io,
    ops::{Index, IndexMut},
};

use noodles_core::position::SequenceIndex;

/// An alignment record quality scores buffer.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct QualityScores(Vec<u8>);

impl QualityScores {
    /// Returns whether there are any scores.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of scores.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns an iterator over the scores.
    pub fn iter(&self) -> impl Iterator<Item = u8> + '_ {
        self.0.iter().copied()
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

impl Extend<u8> for QualityScores {
    fn extend<T: IntoIterator<Item = u8>>(&mut self, iter: T) {
        self.0.extend(iter);
    }
}

impl FromIterator<u8> for QualityScores {
    fn from_iter<T: IntoIterator<Item = u8>>(iter: T) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl From<QualityScores> for Vec<u8> {
    fn from(quality_scores: QualityScores) -> Self {
        quality_scores.0
    }
}

impl crate::alignment::record::QualityScores for &QualityScores {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u8>> + '_> {
        Box::new(self.0.iter().copied().map(Ok))
    }
}
