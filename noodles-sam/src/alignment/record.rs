//! Alignment record.

pub mod sequence;

pub use self::sequence::Sequence;

use self::sequence::Base;

/// An alignment record sequence.
pub trait AlignmentSequence {
    /// Returns the number of bases in the sequence.
    fn len(&self) -> usize;

    /// Returns whether the sequence is empty.
    fn is_empty(&self) -> bool;

    /// Removes all bases from the sequence.
    fn clear(&mut self);

    /// Returns an iterator over the bases in the sequence.
    fn bases(&self) -> Box<dyn Iterator<Item = Base> + '_>;
}
