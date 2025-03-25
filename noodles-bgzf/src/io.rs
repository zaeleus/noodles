//! BGZF I/O.

mod buf_read;
pub mod indexed_reader;
mod read;
pub mod reader;
mod seek;
pub mod writer;

pub use self::{
    buf_read::BufRead, indexed_reader::IndexedReader, read::Read, reader::Reader, seek::Seek,
    writer::Writer,
};
