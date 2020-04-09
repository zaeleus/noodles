pub mod cigar;
pub mod data;
mod flags;
mod header;
mod reader;
mod record;

pub use self::{
    cigar::Cigar, data::Data, flags::Flags, header::Header, reader::Reader, record::Record,
};
