//! Alignment record mapping quality.

use std::{error, fmt, io};

// § 1.4.5 "_MAPQ_" (2023): "A value 255 indicates that the mapping quality is not available."
const MISSING: u8 = 255;

/// An alignment record mapping quality.
///
/// Mapping quality ranges from 0 to 254 (inclusive), where higher is better.
///
/// The value 255 is reserved as a marker for a missing mapping quality.
#[derive(Clone, Copy, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub struct MappingQuality(u8);

impl MappingQuality {
    /// The minimum mapping quality (0).
    pub const MIN: Self = Self(0);

    /// The maximum mapping quality (254).
    pub const MAX: Self = Self(254);

    /// Creates a mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record_buf::MappingQuality;
    /// assert!(MappingQuality::new(8).is_some());
    /// assert!(MappingQuality::new(255).is_none());
    /// ```
    pub const fn new(n: u8) -> Option<Self> {
        if n == MISSING {
            None
        } else {
            Some(Self(n))
        }
    }

    /// Returns the inner value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record_buf::MappingQuality;
    /// let mapping_quality = MappingQuality::new(8).unwrap();
    /// assert_eq!(mapping_quality.get(), 8);
    /// ```
    pub const fn get(&self) -> u8 {
        self.0
    }
}

impl crate::alignment::record::fields::MappingQuality for MappingQuality {
    fn try_to_u8(&self) -> io::Result<u8> {
        Ok(self.get())
    }
}

/// An error returned when a raw alignment record mapping quality fails to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromIntError {
    /// The value is missing.
    Missing,
}

impl error::Error for TryFromIntError {}

impl fmt::Display for TryFromIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing => write!(f, "missing value: {MISSING}"),
        }
    }
}

impl TryFrom<u8> for MappingQuality {
    type Error = TryFromIntError;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        Self::new(n).ok_or(TryFromIntError::Missing)
    }
}

impl From<MappingQuality> for u8 {
    fn from(mapping_quality: MappingQuality) -> Self {
        mapping_quality.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_u8_for_mapping_quality() {
        assert_eq!(MappingQuality::try_from(0), Ok(MappingQuality(0)));
        assert_eq!(MappingQuality::try_from(8), Ok(MappingQuality(8)));
        assert_eq!(MappingQuality::try_from(13), Ok(MappingQuality(13)));
        assert_eq!(MappingQuality::try_from(144), Ok(MappingQuality(144)));
        assert_eq!(MappingQuality::try_from(255), Err(TryFromIntError::Missing));
    }

    #[test]
    fn test_from_mapping_quality_for_u8() {
        assert_eq!(u8::from(MappingQuality(0)), 0);
        assert_eq!(u8::from(MappingQuality(8)), 8);
        assert_eq!(u8::from(MappingQuality(13)), 13);
        assert_eq!(u8::from(MappingQuality(144)), 144);
    }
}
