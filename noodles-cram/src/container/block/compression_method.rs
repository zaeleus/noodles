use std::{convert::TryFrom, error, fmt};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum CompressionMethod {
    None,
    Gzip,
    Bzip2,
    Lzma,
    Rans,
}

impl Default for CompressionMethod {
    fn default() -> Self {
        Self::None
    }
}

#[derive(Debug, Eq, PartialEq)]
pub struct TryFromByteError(u8);

impl error::Error for TryFromByteError {}

impl fmt::Display for TryFromByteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid compression method: expected 0..=4, got {}",
            self.0
        )
    }
}

impl TryFrom<u8> for CompressionMethod {
    type Error = TryFromByteError;

    fn try_from(b: u8) -> Result<Self, Self::Error> {
        match b {
            0 => Ok(Self::None),
            1 => Ok(Self::Gzip),
            2 => Ok(Self::Bzip2),
            3 => Ok(Self::Lzma),
            4 => Ok(Self::Rans),
            _ => Err(TryFromByteError(b)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(CompressionMethod::default(), CompressionMethod::None);
    }

    #[test]
    fn test_try_from() {
        assert_eq!(CompressionMethod::try_from(0), Ok(CompressionMethod::None));
        assert_eq!(CompressionMethod::try_from(1), Ok(CompressionMethod::Gzip));
        assert_eq!(CompressionMethod::try_from(2), Ok(CompressionMethod::Bzip2));
        assert_eq!(CompressionMethod::try_from(3), Ok(CompressionMethod::Lzma));
        assert_eq!(CompressionMethod::try_from(4), Ok(CompressionMethod::Rans));
        assert_eq!(CompressionMethod::try_from(5), Err(TryFromByteError(5)));
    }
}
