//! I/O for variant formats.

mod format;
pub mod reader;
pub mod writer;

pub use self::{
    format::{Compression, Format},
    reader::Reader,
    writer::Writer,
};
