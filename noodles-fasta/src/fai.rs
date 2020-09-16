//! FASTA index (FAI) and fields.

mod reader;
mod record;
mod writer;
mod builder;

pub use self::{reader::Reader, record::Record, writer::Writer, builder::IndexBuilder};
