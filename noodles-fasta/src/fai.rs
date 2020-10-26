//! FASTA index (FAI) and fields.

pub(crate) mod indexer;
mod reader;
mod record;
mod writer;

pub use self::{reader::Reader, record::Record, writer::Writer};
