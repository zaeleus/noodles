//! CRAM container and fields.

pub(crate) mod block;
pub mod block_content_encoder_map;
pub mod compression_header;
mod header;
mod reference_sequence_context;
pub(crate) mod slice;

pub(crate) use self::{
    block::Block, header::Header, reference_sequence_context::ReferenceSequenceContext,
};
pub use self::{
    block_content_encoder_map::BlockContentEncoderMap, compression_header::CompressionHeader,
    slice::Slice,
};

/// A CRAM container.
#[deprecated(since = "0.78.0", note = "Use `cram::io::reader::Container` instead.")]
pub type Container = crate::io::reader::Container;

/// A CRAM container.
#[deprecated(since = "0.78.0", note = "Use `cram::io::reader::Container` instead.")]
pub type DataContainer = crate::io::reader::Container;
