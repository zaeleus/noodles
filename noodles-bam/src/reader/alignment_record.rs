mod cigar;
mod data;
mod mapping_quality;
pub mod quality_scores;
mod read_name;
pub mod sequence;

pub use self::{
    cigar::get_cigar, data::get_data, mapping_quality::get_mapping_quality,
    quality_scores::get_quality_scores, read_name::get_read_name, sequence::get_sequence,
};
