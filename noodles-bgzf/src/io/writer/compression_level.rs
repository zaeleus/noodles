use std::{error, fmt};

const MIN: u8 = 0;

#[cfg(feature = "libdeflate")]
const MAX: u8 = 12;
#[cfg(not(feature = "libdeflate"))]
const MAX: u8 = 9;

/// A DEFLATE compression level.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct CompressionLevel(u8);

impl CompressionLevel {
    /// No compression.
    pub const NONE: Self = Self(0);

    /// A compression level optimized for speed.
    pub const FAST: Self = Self(1);

    /// A compression level optimized for compression rate.
    pub const BEST: Self = Self(MAX);

    /// Creates a compression level.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::io::writer::CompressionLevel;
    /// assert_eq!(CompressionLevel::new(0), Some(CompressionLevel::NONE));
    /// assert!(CompressionLevel::new(255).is_none());
    /// ```
    pub const fn new(n: u8) -> Option<Self> {
        if n <= MAX { Some(Self(n)) } else { None }
    }

    /// Returns the inner value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::io::writer::CompressionLevel;
    /// assert_eq!(CompressionLevel::NONE.get(), 0);
    /// ```
    pub const fn get(&self) -> u8 {
        self.0
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
            Self::Invalid(n) => write!(f, "invalid input: {n}"),
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
    fn test_default() {
        assert_eq!(CompressionLevel::default(), CompressionLevel(6));
    }
}
