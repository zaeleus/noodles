//! I/O for alignment formats.

mod format;
pub mod reader;
pub mod writer;

pub use self::{format::Format, reader::Reader, writer::Writer};
