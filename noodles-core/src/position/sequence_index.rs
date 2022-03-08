use std::ops;

use super::Position;

/// A 1-based index for sequences.
pub trait SequenceIndex<T> {
    /// The output returned.
    type Output: ?Sized;

    /// Returns a reference to the output of the given index.
    fn get(self, sequence: &[T]) -> Option<&Self::Output>;

    /// Returns a mutable reference to the output of the given index.
    fn get_mut(self, sequence: &mut [T]) -> Option<&mut Self::Output>;
}

impl<T> SequenceIndex<T> for Position {
    type Output = T;

    fn get(self, sequence: &[T]) -> Option<&Self::Output> {
        let i = usize::from(self) - 1;
        sequence.get(i)
    }

    fn get_mut(self, sequence: &mut [T]) -> Option<&mut Self::Output> {
        let i = usize::from(self) - 1;
        sequence.get_mut(i)
    }
}

impl<T> SequenceIndex<T> for ops::RangeFrom<Position> {
    type Output = [T];

    fn get(self, sequence: &[T]) -> Option<&Self::Output> {
        let start = usize::from(self.start) - 1;
        sequence.get(start..)
    }

    fn get_mut(self, sequence: &mut [T]) -> Option<&mut Self::Output> {
        let start = usize::from(self.start) - 1;
        sequence.get_mut(start..)
    }
}

impl<T> SequenceIndex<T> for ops::RangeFull {
    type Output = [T];

    fn get(self, sequence: &[T]) -> Option<&Self::Output> {
        Some(sequence)
    }

    fn get_mut(self, sequence: &mut [T]) -> Option<&mut Self::Output> {
        Some(sequence)
    }
}

impl<T> SequenceIndex<T> for ops::RangeInclusive<Position> {
    type Output = [T];

    fn get(self, sequence: &[T]) -> Option<&Self::Output> {
        let start = usize::from(*self.start()) - 1;
        let end = usize::from(*self.end()) - 1;
        sequence.get(start..=end)
    }

    fn get_mut(self, sequence: &mut [T]) -> Option<&mut Self::Output> {
        let start = usize::from(*self.start()) - 1;
        let end = usize::from(*self.end()) - 1;
        sequence.get_mut(start..=end)
    }
}
