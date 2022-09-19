//! I/O for alignment formats.

mod format;
mod reader;
pub mod writer;

pub use self::{format::Format, reader::Reader, writer::Writer};
