//! SAM record and fields.

pub mod data;
pub mod quality_scores;
pub mod sequence;
pub mod template_length;

pub use self::{sequence::Sequence, template_length::TemplateLength};
