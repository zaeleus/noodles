#![warn(missing_docs)]

//! **noodles-bed** handles the reading and writing of the BED (Browser Extensible Data) format.

mod reader;
pub mod record;
mod writer;

pub use self::{reader::Reader, record::Record, writer::Writer};
