//! FASTQ index (FAI) and fields.

mod reader;
mod record;

pub use self::{reader::Reader, record::Record};
