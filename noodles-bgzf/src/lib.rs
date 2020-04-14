pub use self::{block::Block, reader::Reader, writer::Writer};

mod block;
mod gz;
mod reader;
mod writer;

const GZIP_XLEN_SIZE: usize = 2;
const BGZF_XLEN: usize = 6;
pub(crate) const BGZF_HEADER_SIZE: usize = gz::HEADER_SIZE + GZIP_XLEN_SIZE + BGZF_XLEN;
