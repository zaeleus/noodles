//! FAI I/O.

pub(crate) mod reader;
pub(crate) mod writer;

pub use self::{reader::Reader, writer::Writer};
