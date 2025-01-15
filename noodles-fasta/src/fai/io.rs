//! FAI I/O.

mod reader;
pub(crate) mod writer;

pub use self::{reader::Reader, writer::Writer};
