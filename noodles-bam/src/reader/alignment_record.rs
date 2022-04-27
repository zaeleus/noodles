mod cigar;
mod data;
mod quality_scores;
mod read_name;
mod sequence;

pub use self::{
    cigar::get_cigar, data::get_data, quality_scores::get_quality_scores, read_name::get_read_name,
    sequence::get_sequence,
};
