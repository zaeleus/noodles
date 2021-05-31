//! **noodles-csi** handles the reading of the Coordinate-Sorted Index (CSI) format.

mod index;
mod reader;

pub use self::{index::Index, reader::Reader};

static MAGIC_NUMBER: &[u8] = b"CSI\x01";
