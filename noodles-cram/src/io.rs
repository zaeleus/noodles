//! CRAM I/O.

mod bit_reader;
mod bit_writer;
pub mod indexed_reader;
pub mod reader;
pub mod writer;

pub(crate) use self::{bit_reader::BitReader, bit_writer::BitWriter};
pub use self::{indexed_reader::IndexedReader, reader::Reader, writer::Writer};
