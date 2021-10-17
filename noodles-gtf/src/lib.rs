#![warn(missing_docs)]

//! **noodles-gtf** handles the reading and writing of the Gene Transfer Format (GTF).

pub mod line;
mod reader;
pub mod record;
mod writer;

pub use self::{line::Line, reader::Reader, record::Record, writer::Writer};
