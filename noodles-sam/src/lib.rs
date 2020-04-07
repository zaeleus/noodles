pub mod cigar;
mod flags;
mod header;
mod reader;
mod record;

pub use self::{cigar::Cigar, flags::Flags, header::Header, reader::Reader, record::Record};
