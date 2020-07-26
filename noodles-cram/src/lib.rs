mod bit_reader;
mod container;
pub mod crai;
mod data_container;
mod huffman;
mod num;
mod rans;
pub mod reader;
pub mod record;
mod writer;

pub use self::{
    bit_reader::BitReader, container::Container, data_container::DataContainer, reader::Reader,
    record::Record, writer::Writer,
};

static MAGIC_NUMBER: &[u8] = b"CRAM";
