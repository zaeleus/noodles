//! Alignment record field.

mod cigar;
mod data;
mod name;
mod quality_scores;
mod sequence;

pub use self::{
    cigar::Cigar, data::Data, name::Name, quality_scores::QualityScores, sequence::Sequence,
};
