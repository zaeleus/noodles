//! GFF directive value.

pub mod genome_build;
pub mod gff_version;
pub mod sequence_region;

pub use self::{
    genome_build::GenomeBuild, gff_version::GffVersion, sequence_region::SequenceRegion,
};

/// A GFF directive value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Value {
    /// The GFF version (`gff-version`).
    GffVersion(GffVersion),
    /// A reference to a sequence segment (`sequence-region`).
    SequenceRegion(SequenceRegion),
    /// The genome build used for the start and end positions (`genome-build`).
    GenomeBuild(GenomeBuild),
    /// Any other directive value.
    Other(String),
}
