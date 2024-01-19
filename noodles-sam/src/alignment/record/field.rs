//! Alignment record field.

mod cigar;
mod data;
mod flags;
mod mapping_quality;
mod name;
mod position;
mod quality_scores;
mod reference_sequence_id;
mod sequence;
mod template_length;

pub use self::{
    cigar::Cigar, data::Data, flags::Flags, mapping_quality::MappingQuality, name::Name,
    position::Position, quality_scores::QualityScores, reference_sequence_id::ReferenceSequenceId,
    sequence::Sequence, template_length::TemplateLength,
};
