//! Feature record.

pub mod attributes;
mod strand;

use std::io;

use noodles_core::Position;

pub use self::{attributes::Attributes, strand::Strand};
use crate::record::Phase;

/// A feature record.
pub trait Record {
    /// Returns the reference sequence name.
    fn reference_sequence_name(&self) -> &str;

    /// Returns the source.
    fn source(&self) -> &str;

    /// Returns the type.
    fn ty(&self) -> &str;

    /// Returns the feature start.
    fn feature_start(&self) -> io::Result<Position>;

    /// Returns the feature end.
    fn feature_end(&self) -> io::Result<Position>;

    /// Returns the score.
    fn score(&self) -> Option<io::Result<f32>>;

    /// Returns the strand.
    fn strand(&self) -> io::Result<Strand>;

    /// Returns the phase.
    fn phase(&self) -> Option<io::Result<Phase>>;

    /// Returns the attributes.
    fn attributes(&self) -> Box<dyn Attributes + '_>;
}
