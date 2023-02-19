//! CRAM container block compression method.

/// A CRAM container block compression method.
///
/// The compression method is associated with a block to identify how its data is compressed.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[non_exhaustive]
pub enum CompressionMethod {
    /// Uncompressed.
    #[default]
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(CompressionMethod::default(), CompressionMethod::None);
    }
}
