//! I/O for variant formats.

mod format;
pub mod reader;

pub use self::{
    format::{Compression, Format},
    reader::Reader,
};
