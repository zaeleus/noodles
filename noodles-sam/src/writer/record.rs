mod cigar;
mod data;
mod position;
mod quality_scores;
mod sequence;

pub use self::{
    cigar::write_cigar, data::write_data, position::write_position,
    quality_scores::write_quality_scores, sequence::write_sequence,
};

const MISSING: u8 = b'*';
