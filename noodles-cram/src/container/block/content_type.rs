//! CRAM container block content type.

/// A CRAM container block content type.
///
/// The content type is associated with a block to identify its data.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum ContentType {
    /// CRAM file header.
    ///
    /// The CRAM file header contains the SAM header.
    FileHeader,
    /// Compression header.
    ///
    /// There is a single compression header associated with each container. It contains the
    /// preservation map, data series encodings, and tag encoding map.
    CompressionHeader,
    /// Slice header.
    ///
    /// A slice typically has a contiguous region of alignments records.
    SliceHeader,
    /// Reserved.
    Reserved,
    /// External data.
    ///
    /// External data is used in slices. There is at least one external block, typically
    /// corresponding to a byte stream of external data.
    ExternalData,
    /// Core data.
    ///
    /// Core data is used in slices. It is a single block with a bit stream of any non-external
    /// encoding.
    CoreData,
}
