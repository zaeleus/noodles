//! **noodles-csi** handles the reading of the Coordinate-Sorted Index (CSI) format.

mod index;
mod reader;
mod writer;

pub use self::{index::Index, reader::Reader, writer::Writer};

static MAGIC_NUMBER: &[u8] = b"CSI\x01";
