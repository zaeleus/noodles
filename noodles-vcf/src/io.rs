//! VCF I/O.

mod compression_method;
pub mod indexed_reader;
pub mod reader;

pub use self::{
    compression_method::CompressionMethod, indexed_reader::IndexedReader, reader::Reader,
};
