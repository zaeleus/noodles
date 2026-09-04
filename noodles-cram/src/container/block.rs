//! CRAM container block and fields.

mod compression_method;
mod content_type;

pub use self::{compression_method::CompressionMethod, content_type::ContentType};

/// A CRAM container block content ID.
///
/// This associates an external data block with a data series.
pub type ContentId = i32;
