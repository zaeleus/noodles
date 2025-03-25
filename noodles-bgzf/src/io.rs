//! BGZF I/O.

mod buf_read;
pub mod indexed_reader;
mod multithreaded_reader;
mod read;
pub mod reader;
mod seek;
pub mod writer;

pub use self::{
    buf_read::BufRead, indexed_reader::IndexedReader, multithreaded_reader::MultithreadedReader,
    read::Read, reader::Reader, seek::Seek, writer::Writer,
};
