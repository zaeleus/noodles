use std::ops;

use super::Position;

/// A 1-based index for sequences.
pub trait SequenceIndex {
    /// The output returned.
    type Output: ?Sized;

    /// Returns a reference to the output of the given index.
    fn get(self, sequence: &[u8]) -> Option<&Self::Output>;
}

impl SequenceIndex for Position {
    type Output = u8;

    fn get(self, sequence: &[u8]) -> Option<&Self::Output> {
        let i = usize::from(self) - 1;
        sequence.get(i)
    }
}

impl SequenceIndex for ops::RangeFrom<Position> {
    type Output = [u8];

    fn get(self, sequence: &[u8]) -> Option<&Self::Output> {
        let start = usize::from(self.start) - 1;
        sequence.get(start..)
    }
}

impl SequenceIndex for ops::RangeFull {
    type Output = [u8];

    fn get(self, sequence: &[u8]) -> Option<&Self::Output> {
        Some(sequence)
    }
}

impl SequenceIndex for ops::RangeInclusive<Position> {
    type Output = [u8];

    fn get(self, sequence: &[u8]) -> Option<&Self::Output> {
        let start = usize::from(*self.start()) - 1;
        let end = usize::from(*self.end()) - 1;
        sequence.get(start..=end)
    }
}
