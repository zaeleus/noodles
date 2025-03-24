//! BGZF I/O.

mod buf_read;
mod read;
pub mod reader;
mod seek;

pub use self::{buf_read::BufRead, read::Read, reader::Reader, seek::Seek};
