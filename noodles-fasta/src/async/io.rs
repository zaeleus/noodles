//! Async FASTA I/O.

pub(crate) mod reader;
pub mod writer;

pub use self::{reader::Reader, writer::Writer};
