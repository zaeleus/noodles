#![warn(missing_docs)]

//! **noodles-cram** handles the reading and writing of the CRAM format.

#[cfg(feature = "async")]
mod r#async;

mod bit_reader;
mod bit_writer;
mod codecs;
pub(crate) mod container;
pub mod crai;
pub mod data_container;
pub mod file_definition;
mod huffman;
mod indexer;
mod num;
pub mod reader;
pub mod record;
pub mod writer;

pub use self::{
    data_container::DataContainer, file_definition::FileDefinition, indexer::index, reader::Reader,
    record::Record, writer::Writer,
};

#[cfg(feature = "async")]
pub use self::r#async::Reader as AsyncReader;

pub(crate) use self::{bit_reader::BitReader, bit_writer::BitWriter};

static MAGIC_NUMBER: &[u8] = b"CRAM";
