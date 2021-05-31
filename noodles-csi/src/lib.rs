//! **noodles-csi** handles the reading and writing of the coordinate-sorted index (CSI) format.

mod index;
mod reader;
mod writer;

pub use self::{index::Index, reader::Reader, writer::Writer};

static MAGIC_NUMBER: &[u8] = b"CSI\x01";
