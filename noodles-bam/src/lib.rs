pub mod bai;
pub mod cigar;
pub mod data;
mod quality;
pub mod reader;
mod record;
mod reference;
pub mod sequence;

pub use self::{
    cigar::Cigar,
    data::Data,
    quality::Quality,
    reader::{Reader, Records, References},
    record::Record,
    reference::Reference,
    sequence::Sequence,
};

pub static MAGIC_NUMBER: &[u8] = b"BAM\x01";
