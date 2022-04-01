mod cigar;
mod data;
mod quality_scores;
mod sequence;

pub use self::{
    cigar::write_cigar, data::write_data, quality_scores::write_quality_scores,
    sequence::write_sequence,
};

const MISSING: u8 = b'*';
