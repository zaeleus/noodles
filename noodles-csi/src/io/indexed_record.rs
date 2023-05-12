use noodles_core::{region::Interval, Position};

/// An indexed record.
pub trait IndexedRecord {
    /// Returns the reference sequence name.
    fn indexed_reference_sequence_name(&self) -> &str;

    /// Returns the start position.
    fn indexed_start_position(&self) -> Position;

    /// Returns the end position.
    fn indexed_end_position(&self) -> Position;

    /// Returns the start and end positions as an [`Interval`].
    fn indexed_interval(&self) -> Interval {
        (self.indexed_start_position()..=self.indexed_end_position()).into()
    }
}
