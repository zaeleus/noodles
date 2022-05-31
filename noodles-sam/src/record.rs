//! SAM record and fields.

pub mod cigar;
pub mod data;
mod flags;
pub mod mapping_quality;
pub mod quality_scores;
pub mod read_name;
pub mod reference_sequence_name;
pub mod sequence;
pub mod template_length;

pub use self::{
    cigar::Cigar, data::Data, flags::Flags, mapping_quality::MappingQuality,
    quality_scores::QualityScores, read_name::ReadName,
    reference_sequence_name::ReferenceSequenceName, sequence::Sequence,
    template_length::TemplateLength,
};
