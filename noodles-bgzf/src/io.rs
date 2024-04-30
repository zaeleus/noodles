//! BGZF I/O.

mod buf_read;
mod read;

pub use self::{buf_read::BufRead, read::Read};
