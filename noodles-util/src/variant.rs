//! I/O for variant formats.

mod compression_method;
mod format;
pub mod indexed_reader;
pub mod reader;
pub mod writer;

pub use self::{
    compression_method::CompressionMethod, format::Format, indexed_reader::IndexedReader,
    reader::Reader, writer::Writer,
};
