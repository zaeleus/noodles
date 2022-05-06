//! Alignment record.

pub mod cigar;
pub mod data;
mod flags;
pub mod mapping_quality;
pub mod quality_scores;
pub mod read_name;
pub mod sequence;

pub use self::{
    cigar::Cigar, data::Data, flags::Flags, mapping_quality::MappingQuality,
    quality_scores::QualityScores, read_name::ReadName, sequence::Sequence,
};
