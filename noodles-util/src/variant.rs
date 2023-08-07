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

/// A variant compression method.
#[deprecated(since = "0.20.0", note = "Use `CompressionMethod` instead.")]
pub type Compression = CompressionMethod;
