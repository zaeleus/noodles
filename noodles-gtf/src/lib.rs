#![warn(missing_docs)]

//! **noodles-gtf** handles the reading of the Gene Transfer Format (GTF).

mod reader;
pub mod record;

pub use self::{reader::Reader, record::Record};
