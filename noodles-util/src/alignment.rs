//! I/O for alignment formats.

mod format;
pub mod indexed_reader;
pub mod reader;
pub mod writer;

pub use self::{format::Format, indexed_reader::IndexedReader, reader::Reader, writer::Writer};
