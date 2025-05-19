//! FASTQ I/O.

mod indexer;
pub mod reader;
pub mod writer;

pub use self::{indexer::Indexer, reader::Reader, writer::Writer};
