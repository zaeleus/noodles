pub mod compression_header;
pub mod slice;

pub use self::{compression_header::read_compression_header, slice::read_slice};
