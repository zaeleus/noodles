//! BGZF I/O.

mod block;
mod buf_read;
pub mod indexed_reader;
mod multithreaded_reader;
pub mod multithreaded_writer;
mod read;
pub mod reader;
mod seek;
pub mod writer;

pub(crate) use self::block::Block;
pub use self::{
    buf_read::BufRead, indexed_reader::IndexedReader, multithreaded_reader::MultithreadedReader,
    multithreaded_writer::MultithreadedWriter, read::Read, reader::Reader, seek::Seek,
    writer::Writer,
};
