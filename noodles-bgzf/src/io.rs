//! BGZF I/O.

mod buf_read;
mod read;
mod seek;

pub use self::{buf_read::BufRead, read::Read, seek::Seek};
