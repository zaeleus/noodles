//! SAM record and fields.

pub mod data;
mod flags;
pub mod mapping_quality;
pub mod name;
pub mod quality_scores;
pub mod reference_sequence_name;
pub mod sequence;
pub mod template_length;

pub use self::{
    flags::Flags, mapping_quality::MappingQuality, name::Name, quality_scores::QualityScores,
    reference_sequence_name::ReferenceSequenceName, sequence::Sequence,
    template_length::TemplateLength,
};
