#![warn(missing_docs)]

//! **noodles-gtf** handles the reading of the Gene Transfer Format (GTF).

pub mod line;
mod reader;
pub mod record;

pub use self::{line::Line, reader::Reader, record::Record};
