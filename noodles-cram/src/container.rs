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
#[deprecated(since = "0.78.0", note = "Use `cram::Container` instead.")]
pub type DataContainer = Container;

/// A CRAM container.
pub struct Container {
    compression_header: CompressionHeader,
    slices: Vec<Slice>,
}

impl Container {
    pub(crate) fn new(compression_header: CompressionHeader, slices: Vec<Slice>) -> Self {
        Self {
            compression_header,
            slices,
        }
    }

    /// Returns the compression header.
    pub fn compression_header(&self) -> &CompressionHeader {
        &self.compression_header
    }

    /// Returns a list of slices.
    pub fn slices(&self) -> &[Slice] {
        &self.slices
    }
}
