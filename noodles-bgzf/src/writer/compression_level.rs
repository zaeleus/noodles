use std::{error, fmt};

const MIN: u8 = 0;

#[cfg(feature = "libdeflate")]
const MAX: u8 = 12;
#[cfg(not(feature = "libdeflate"))]
const MAX: u8 = 9;

/// A DEFLATE compression level.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CompressionLevel(u8);

impl CompressionLevel {
    /// Returns a compression level to disable compression.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::writer::CompressionLevel;
    /// let compression_level = CompressionLevel::none();
    /// ```
    pub fn none() -> Self {
        Self(0)
    }

    /// Returns a compression level optimized for speed.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::writer::CompressionLevel;
    /// let compression_level = CompressionLevel::fast();
    /// ```
    pub fn fast() -> Self {
        Self(1)
    }

    /// Returns a compression level optimized for compression rate.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::writer::CompressionLevel;
    /// let compression_level = CompressionLevel::best();
    /// ```
    pub fn best() -> Self {
        Self(MAX)
    }
}

impl Default for CompressionLevel {
    fn default() -> Self {
        Self(6)
    }
}

/// An error returned when a raw DEFLATE compression level fails to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromU8Error {
    Invalid(u8),
}

impl error::Error for TryFromU8Error {}

impl fmt::Display for TryFromU8Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(n) => write!(f, "invalid input: {}", n),
        }
    }
}

impl TryFrom<u8> for CompressionLevel {
    type Error = TryFromU8Error;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        match n {
            MIN..=MAX => Ok(Self(n)),
            _ => Err(TryFromU8Error::Invalid(n)),
        }
    }
}

impl From<CompressionLevel> for u8 {
    fn from(compression_level: CompressionLevel) -> Self {
        compression_level.0
    }
}

#[cfg(feature = "libdeflate")]
impl From<CompressionLevel> for libdeflater::CompressionLvl {
    fn from(compression_level: CompressionLevel) -> Self {
        // SAFETY: The raw value is guaranteed to be between MIN and MAX, inclusive.
        Self::new(i32::from(u8::from(compression_level))).unwrap()
    }
}

#[cfg(not(feature = "libdeflate"))]
impl From<CompressionLevel> for flate2::Compression {
    fn from(compression_level: CompressionLevel) -> Self {
        Self::new(u32::from(u8::from(compression_level)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_none() {
        assert_eq!(CompressionLevel::none(), CompressionLevel(0));
    }

    #[test]
    fn test_fast() {
        assert_eq!(CompressionLevel::fast(), CompressionLevel(1));
    }

    #[test]
    fn test_best() {
        assert_eq!(CompressionLevel::best(), CompressionLevel(MAX));
    }

    #[test]
    fn test_default() {
        assert_eq!(CompressionLevel::default(), CompressionLevel(6));
    }
}
