//! Variant record.

mod alternate_bases;
mod filters;
mod ids;
pub mod info;
pub mod samples;

use std::io;

use noodles_core::Position;

pub use self::{
    alternate_bases::AlternateBases, filters::Filters, ids::Ids, info::Info, samples::Samples,
};
use crate::Header;

/// A variant record.
pub trait Record {
    /// Returns the reference sequence name.
    fn reference_sequence_name(&self, header: &Header) -> &str;

    /// Returns the position.
    fn position(&self) -> Option<io::Result<Position>>;

    /// Returns the IDs.
    fn ids(&self) -> Box<dyn Ids + '_>;

    /// Returns the reference bases.
    fn reference_bases(&self) -> &str;

    /// Returns the alternate bases.
    fn alternate_bases(&self) -> Box<dyn AlternateBases + '_>;

    /// Returns the quality scores.
    fn quality_score(&self) -> Option<io::Result<f32>>;

    /// Returns the filters.
    fn filters(&self) -> Box<dyn Filters + '_>;

    /// Return the info fields.
    fn info(&self) -> Box<dyn Info + '_>;

    /// Returns the samples.
    fn samples(&self) -> Box<dyn Samples + '_>;
}
