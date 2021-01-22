//! FASTA index (FAI) and fields.

mod reader;
mod record;
mod writer;

pub use self::{reader::Reader, record::Record, writer::Writer};

/// A FASTA index.
pub type Index = Vec<Record>;
