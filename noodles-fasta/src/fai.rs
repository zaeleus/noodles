//! FASTA index (FAI) and fields.

mod builder;
mod reader;
mod record;
mod writer;

pub use self::{builder::IndexBuilder, reader::Reader, record::Record, writer::Writer};
