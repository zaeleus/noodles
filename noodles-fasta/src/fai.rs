//! FASTA index (FAI) and fields.

mod indexer;
mod reader;
mod record;
mod writer;

pub use self::{indexer::Indexer, reader::Reader, record::Record, writer::Writer};
