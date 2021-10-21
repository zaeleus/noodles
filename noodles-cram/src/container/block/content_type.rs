//! CRAM container block content type.

use std::{error, fmt};

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

/// An error returned when a raw content type fails to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromByteError(u8);

impl error::Error for TryFromByteError {}

impl fmt::Display for TryFromByteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid content type: expected 0..=5, got {}", self.0)
    }
}

impl TryFrom<u8> for ContentType {
    type Error = TryFromByteError;

    fn try_from(b: u8) -> Result<Self, Self::Error> {
        match b {
            0 => Ok(Self::FileHeader),
            1 => Ok(Self::CompressionHeader),
            2 => Ok(Self::SliceHeader),
            3 => Ok(Self::Reserved),
            4 => Ok(Self::ExternalData),
            5 => Ok(Self::CoreData),
            _ => Err(TryFromByteError(b)),
        }
    }
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
    fn test_try_from_u8_for_content_type() {
        assert_eq!(ContentType::try_from(0), Ok(ContentType::FileHeader));
        assert_eq!(ContentType::try_from(1), Ok(ContentType::CompressionHeader));
        assert_eq!(ContentType::try_from(2), Ok(ContentType::SliceHeader));
        assert_eq!(ContentType::try_from(3), Ok(ContentType::Reserved));
        assert_eq!(ContentType::try_from(4), Ok(ContentType::ExternalData));
        assert_eq!(ContentType::try_from(5), Ok(ContentType::CoreData));
        assert_eq!(ContentType::try_from(6), Err(TryFromByteError(6)));
    }

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
