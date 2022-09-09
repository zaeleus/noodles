//! CRAM container block compression method.

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
