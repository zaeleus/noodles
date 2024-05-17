//! GFF I/O.

pub(crate) mod reader;
mod writer;

pub use self::{reader::Reader, writer::Writer};
