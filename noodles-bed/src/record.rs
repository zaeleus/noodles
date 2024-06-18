mod fields;

use std::io;

use noodles_core::Position;

use self::fields::Fields;

/// A BED record.
#[derive(Clone, Default, Eq, PartialEq)]
pub struct Record(Fields);

impl Record {
    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &[u8] {
        self.0.reference_sequence_name()
    }

    /// Returns the feature start.
    pub fn feature_start(&self) -> io::Result<Position> {
        self.0.feature_start()
    }

    /// Returns the feature end.
    pub fn feature_end(&self) -> io::Result<Position> {
        self.0.feature_end()
    }

    /// Returns the name.
    pub fn name(&self) -> Option<&[u8]> {
        self.0.name()
    }

    /// Returns the score.
    pub fn score(&self) -> Option<&[u8]> {
        self.0.score()
    }

    /// Returns the strand.
    pub fn strand(&self) -> Option<&[u8]> {
        self.0.strand()
    }
}
