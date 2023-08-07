//! I/O for variant formats.

mod compression;
mod format;
pub mod indexed_reader;
pub mod reader;
pub mod writer;

pub use self::{
    compression::Compression, format::Format, indexed_reader::IndexedReader, reader::Reader,
    writer::Writer,
};
