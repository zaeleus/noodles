//! BAM I/O.

pub mod indexed_reader;
pub mod reader;

pub use self::{indexed_reader::IndexedReader, reader::Reader};
