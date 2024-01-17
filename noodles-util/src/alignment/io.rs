//! Alignment format I/O.

mod compression_method;
mod format;
pub mod indexed_reader;
pub mod reader;

pub use self::{
    compression_method::CompressionMethod, format::Format, indexed_reader::IndexedReader,
    reader::Reader,
};
