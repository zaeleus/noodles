//! Alignment record.

pub mod quality_scores;
pub mod sequence;

pub use self::{quality_scores::QualityScores, sequence::Sequence};

use std::io;

use self::{quality_scores::Score, sequence::Base};

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

/// Alignment record quality scores.
pub trait AlignmentQualityScores {
    /// Returns the number of scores.
    fn len(&self) -> usize;

    /// Returns whether there are any scores.
    fn is_empty(&self) -> bool;

    /// Removes all scores.
    fn clear(&mut self);

    /// Returns an iterator over the scores.
    fn scores(&self) -> Box<dyn Iterator<Item = io::Result<Score>> + '_>;
}
