mod bit_reader;
mod container;
pub mod crai;
mod huffman;
mod num;
mod rans;
mod reader;
pub mod record;
mod writer;

pub use self::{
    bit_reader::BitReader, container::Container, reader::Reader, record::Record, writer::Writer,
};

static MAGIC_NUMBER: &[u8] = b"CRAM";
