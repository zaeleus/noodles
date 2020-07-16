use std::{convert::TryFrom, error, fmt};

#[derive(Debug, Eq, PartialEq)]
pub struct TryFromByteError(u8);

impl fmt::Display for TryFromByteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid content type: expected 0..=5, got {}", self.0)
    }
}

impl error::Error for TryFromByteError {}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum ContentType {
    FileHeader,
    CompressionHeader,
    SliceHeader,
    Reserved,
    ExternalData,
    CoreData,
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

#[cfg(test)]
mod tests {
    use std::convert::TryFrom;

    use super::*;

    #[test]
    fn test_try_from() {
        assert_eq!(ContentType::try_from(0), Ok(ContentType::FileHeader));
        assert_eq!(ContentType::try_from(1), Ok(ContentType::CompressionHeader));
        assert_eq!(ContentType::try_from(2), Ok(ContentType::SliceHeader));
        assert_eq!(ContentType::try_from(3), Ok(ContentType::Reserved));
        assert_eq!(ContentType::try_from(4), Ok(ContentType::ExternalData));
        assert_eq!(ContentType::try_from(5), Ok(ContentType::CoreData));
        assert_eq!(ContentType::try_from(6), Err(TryFromByteError(6)));
    }
}
