//! SAM record and fields.

pub mod data;
pub mod quality_scores;
pub mod sequence;
pub mod template_length;

pub use self::{
    quality_scores::QualityScores, sequence::Sequence, template_length::TemplateLength,
};
