pub mod cigar;
mod flags;
mod header;
mod reader;

pub use self::{cigar::Cigar, flags::Flags, header::Header, reader::Reader};
