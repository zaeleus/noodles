//! **noodles-bcf** handles the reading and writing of the BCF format.

pub mod header;
mod reader;
mod record;

pub use self::{reader::Reader, record::Record};
