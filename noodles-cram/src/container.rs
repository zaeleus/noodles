//! CRAM container and fields.

pub(crate) mod block;
pub mod block_content_encoder_map;
pub mod compression_header;
mod header;
mod reference_sequence_context;
pub(crate) mod slice;

pub use self::{
    block_content_encoder_map::BlockContentEncoderMap, compression_header::CompressionHeader,
};
pub(crate) use self::{
    header::Header,
    reference_sequence_context::{
        Context as ReferenceSequenceContextInner, ReferenceSequenceContext,
    },
};
