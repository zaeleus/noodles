//! SAM I/O.

mod compression_method;
pub mod indexed_reader;
pub mod reader;
pub mod writer;

pub use self::{
    compression_method::CompressionMethod, indexed_reader::IndexedReader, reader::Reader,
    writer::Writer,
};
