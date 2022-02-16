//! CRAM container block compression method.

use std::{error, fmt};

/// A CRAM container block compression method.
///
/// The compression method is associated with a block to identify how its data is compressed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum CompressionMethod {
    /// Uncompressed.
    None,
    /// gzip.
    Gzip,
    /// bzip2.
    Bzip2,
    /// Lempel-Ziv-Markov chain algorithm (LZMA).
    Lzma,
    /// Ranged asymmetric numeral systems (4 states, 8-bit renormalization; rANS 4x8).
    Rans4x8,
    /// Ranged asymmetric numeral systems (N states, 16-bit renormalization; rANS Nx16).
    RansNx16,
    /// Adaptive arithmetic coding.
    AdaptiveArithmeticCoding,
    /// fqzcomp.
    Fqzcomp,
    /// Name tokenization codec.
    NameTokenizer,
}

impl Default for CompressionMethod {
    fn default() -> Self {
        Self::None
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromByteError(u8);

impl error::Error for TryFromByteError {}

impl fmt::Display for TryFromByteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid compression method: expected 0..=8, got {}",
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
            4 => Ok(Self::Rans4x8),
            5 => Ok(Self::RansNx16),
            6 => Ok(Self::AdaptiveArithmeticCoding),
            7 => Ok(Self::Fqzcomp),
            8 => Ok(Self::NameTokenizer),
            _ => Err(TryFromByteError(b)),
        }
    }
}

impl From<CompressionMethod> for u8 {
    fn from(compression_method: CompressionMethod) -> Self {
        match compression_method {
            CompressionMethod::None => 0,
            CompressionMethod::Gzip => 1,
            CompressionMethod::Bzip2 => 2,
            CompressionMethod::Lzma => 3,
            CompressionMethod::Rans4x8 => 4,
            CompressionMethod::RansNx16 => 5,
            CompressionMethod::AdaptiveArithmeticCoding => 6,
            CompressionMethod::Fqzcomp => 7,
            CompressionMethod::NameTokenizer => 8,
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
    fn test_try_from_u8_for_compression_method() {
        assert_eq!(CompressionMethod::try_from(0), Ok(CompressionMethod::None));
        assert_eq!(CompressionMethod::try_from(1), Ok(CompressionMethod::Gzip));
        assert_eq!(CompressionMethod::try_from(2), Ok(CompressionMethod::Bzip2));
        assert_eq!(CompressionMethod::try_from(3), Ok(CompressionMethod::Lzma));
        assert_eq!(
            CompressionMethod::try_from(4),
            Ok(CompressionMethod::Rans4x8)
        );
        assert_eq!(
            CompressionMethod::try_from(5),
            Ok(CompressionMethod::RansNx16)
        );
        assert_eq!(
            CompressionMethod::try_from(6),
            Ok(CompressionMethod::AdaptiveArithmeticCoding)
        );
        assert_eq!(
            CompressionMethod::try_from(7),
            Ok(CompressionMethod::Fqzcomp)
        );
        assert_eq!(
            CompressionMethod::try_from(8),
            Ok(CompressionMethod::NameTokenizer)
        );
        assert_eq!(CompressionMethod::try_from(9), Err(TryFromByteError(9)));
    }

    #[test]
    fn test_from_compression_method_for_u8() {
        assert_eq!(u8::from(CompressionMethod::None), 0);
        assert_eq!(u8::from(CompressionMethod::Gzip), 1);
        assert_eq!(u8::from(CompressionMethod::Bzip2), 2);
        assert_eq!(u8::from(CompressionMethod::Lzma), 3);
        assert_eq!(u8::from(CompressionMethod::Rans4x8), 4);
        assert_eq!(u8::from(CompressionMethod::RansNx16), 5);
        assert_eq!(u8::from(CompressionMethod::AdaptiveArithmeticCoding), 6);
        assert_eq!(u8::from(CompressionMethod::Fqzcomp), 7);
        assert_eq!(u8::from(CompressionMethod::NameTokenizer), 8);
    }
}
