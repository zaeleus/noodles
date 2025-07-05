//! BAM I/O.

pub mod indexed_reader;
pub mod reader;
pub mod writer;

pub use self::{indexed_reader::IndexedReader, reader::Reader, writer::Writer};

pub(crate) const MAGIC_NUMBER: [u8; 4] = *b"BAM\x01";
