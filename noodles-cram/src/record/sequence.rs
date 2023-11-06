use std::ops::{Index, IndexMut};

use noodles_core::position::SequenceIndex;

/// RAW CRAM record sequence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Sequence(Vec<u8>);

impl Sequence {
    /// Returns whether there are any bases.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of bases.
    pub fn len(&self) -> usize {
        self.0.len()
    }
}

impl AsRef<[u8]> for Sequence {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<Vec<u8>> for Sequence {
    fn as_mut(&mut self) -> &mut Vec<u8> {
        &mut self.0
    }
}

impl From<Vec<u8>> for Sequence {
    fn from(bases: Vec<u8>) -> Self {
        Self(bases)
    }
}

impl<I> Index<I> for Sequence
where
    I: SequenceIndex<u8>,
{
    type Output = I::Output;

    fn index(&self, index: I) -> &Self::Output {
        index.index(&self.0)
    }
}

impl<I> IndexMut<I> for Sequence
where
    I: SequenceIndex<u8>,
{
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        index.index_mut(&mut self.0)
    }
}

impl From<Sequence> for Vec<u8> {
    fn from(sequence: Sequence) -> Self {
        sequence.0
    }
}
