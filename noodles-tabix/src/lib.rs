pub mod index;
mod reader;
mod writer;

pub use self::{index::Index, reader::Reader, writer::Writer};

static MAGIC_NUMBER: &[u8] = b"TBI\x01";
