pub mod index;
mod reader;

pub use self::{index::Index, reader::Reader};

static MAGIC_NUMBER: &[u8] = b"TBI\x01";
