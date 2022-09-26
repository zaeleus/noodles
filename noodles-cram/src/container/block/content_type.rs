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
    /// There is a single compression header associated with each data container. It contains the
    /// preservation map, data series encoding map, and tag encoding map.
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

impl From<ContentType> for u8 {
    fn from(content_type: ContentType) -> Self {
        match content_type {
            ContentType::FileHeader => 0,
            ContentType::CompressionHeader => 1,
            ContentType::SliceHeader => 2,
            ContentType::Reserved => 3,
            ContentType::ExternalData => 4,
            ContentType::CoreData => 5,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_content_type_for_u8() {
        assert_eq!(u8::from(ContentType::FileHeader), 0);
        assert_eq!(u8::from(ContentType::CompressionHeader), 1);
        assert_eq!(u8::from(ContentType::SliceHeader), 2);
        assert_eq!(u8::from(ContentType::Reserved), 3);
        assert_eq!(u8::from(ContentType::ExternalData), 4);
        assert_eq!(u8::from(ContentType::CoreData), 5);
    }
}
