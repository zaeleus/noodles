//! **noodles-bcf** handles the reading and writing of the BCF format.

pub mod header;
mod reader;
mod record;
mod writer;

pub use self::{reader::Reader, record::Record, writer::Writer};

static MAGIC_NUMBER: &[u8] = b"BCF";
