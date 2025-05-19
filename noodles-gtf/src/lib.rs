//! **noodles-gtf** handles the reading and writing of the Gene Transfer Format (GTF).

pub mod io;
pub mod line;
mod line_buf;
pub mod record;

pub use self::{line::Line, line_buf::LineBuf, record::Record};
