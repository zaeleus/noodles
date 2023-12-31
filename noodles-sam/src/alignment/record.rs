//! Alignment record.

mod cigar;
pub mod data;
mod flags;
mod mapping_quality;
mod position;
mod quality_scores;
mod read_name;
mod reference_sequence_id;
mod sequence;
mod template_length;

pub use self::{
    cigar::Cigar, data::Data, flags::Flags, mapping_quality::MappingQuality, position::Position,
    quality_scores::QualityScores, read_name::ReadName, reference_sequence_id::ReferenceSequenceId,
    sequence::Sequence, template_length::TemplateLength,
};
