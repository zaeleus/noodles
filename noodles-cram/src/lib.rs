mod bit_reader;
mod bit_writer;
pub mod container;
pub mod crai;
mod data_container;
pub mod file_definition;
mod huffman;
mod num;
mod rans;
pub mod reader;
pub mod record;
pub mod writer;

pub use self::{
    bit_reader::BitReader, bit_writer::BitWriter, container::Container,
    data_container::DataContainer, file_definition::FileDefinition, reader::Reader, record::Record,
    writer::Writer,
};

static MAGIC_NUMBER: &[u8] = b"CRAM";
