//! Feature record.

mod other_fields;

use std::io;

use noodles_core::Position;

pub use self::other_fields::OtherFields;

/// A feature record.
pub trait Record {
    /// Return the number of standard fields.
    fn standard_field_count(&self) -> usize;

    /// Returns the reference sequence name.
    fn reference_sequence_name(&self) -> &[u8];

    /// Returns the feature start.
    fn feature_start(&self) -> io::Result<Position>;

    /// Returns the feature end.
    fn feature_end(&self) -> Option<io::Result<Position>>;
}
