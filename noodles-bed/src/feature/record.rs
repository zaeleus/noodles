//! Feature record.

pub mod other_fields;
mod strand;

use std::io;

use bstr::BStr;
use noodles_core::Position;

pub use self::{other_fields::OtherFields, strand::Strand};

/// A feature record.
pub trait Record<const N: usize> {
    /// Return the number of standard fields.
    fn standard_field_count(&self) -> usize {
        N
    }

    /// Returns the reference sequence name.
    fn reference_sequence_name(&self) -> &BStr;

    /// Returns the feature start.
    fn feature_start(&self) -> io::Result<Position>;

    /// Returns the feature end.
    fn feature_end(&self) -> Option<io::Result<Position>>;

    /// Returns the name.
    fn name(&self) -> Option<Option<&BStr>>;

    /// Returns the score.
    fn score(&self) -> Option<io::Result<u16>>;

    /// Returns the strand.
    fn strand(&self) -> Option<io::Result<Option<Strand>>>;

    /// Returns the other fields.
    fn other_fields(&self) -> Box<dyn OtherFields + '_>;
}
