pub mod bai;
pub mod cigar;
pub mod data;
mod quality;
pub mod reader;
mod record;
mod reference;
pub mod sequence;
mod writer;

pub use self::{
    cigar::Cigar,
    data::Data,
    quality::Quality,
    reader::{Reader, Records},
    record::Record,
    reference::Reference,
    sequence::Sequence,
    writer::Writer,
};

pub static MAGIC_NUMBER: &[u8] = b"BAM\x01";
