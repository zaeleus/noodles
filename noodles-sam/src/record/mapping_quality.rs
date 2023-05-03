//! SAM record mapping quality.

use std::{error, fmt, num, str::FromStr};

/// The raw value of a missing mapping quality.
pub const MISSING: u8 = 255;

/// A SAM record mapping quality.
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

    /// Creates a mapping quality if the given value is not missing.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::MappingQuality;
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
    /// use noodles_sam::record::MappingQuality;
    /// let mapping_quality = MappingQuality::new(8).unwrap();
    /// assert_eq!(mapping_quality.get(), 8);
    /// ```
    pub const fn get(&self) -> u8 {
        self.0
    }
}

impl fmt::Display for MappingQuality {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// An error returned when a raw SAM record mapping quality fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input failed to parse as an integer.
    Parse(num::ParseIntError),
    /// The value is missing.
    Missing,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Parse(e) => Some(e),
            Self::Missing => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(_) => f.write_str("parse error"),
            Self::Missing => write!(f, "missing value: {MISSING}"),
        }
    }
}

impl FromStr for MappingQuality {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let n: u8 = s.parse().map_err(ParseError::Parse)?;
        Self::try_from(n)
    }
}

impl TryFrom<u8> for MappingQuality {
    type Error = ParseError;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        Self::new(n).ok_or(ParseError::Missing)
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
    fn test_fmt() {
        assert_eq!(MappingQuality(0).to_string(), "0");
        assert_eq!(MappingQuality(8).to_string(), "8");
        assert_eq!(MappingQuality(13).to_string(), "13");
        assert_eq!(MappingQuality(144).to_string(), "144");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("0".parse(), Ok(MappingQuality(0)));
        assert_eq!("8".parse(), Ok(MappingQuality(8)));
        assert_eq!("13".parse(), Ok(MappingQuality(13)));
        assert_eq!("144".parse(), Ok(MappingQuality(144)));

        assert!(matches!(
            "".parse::<MappingQuality>(),
            Err(ParseError::Parse(_))
        ));
        assert!(matches!(
            "256".parse::<MappingQuality>(),
            Err(ParseError::Parse(_))
        ));
        assert_eq!("255".parse::<MappingQuality>(), Err(ParseError::Missing));
    }

    #[test]
    fn test_try_from_u8_for_mapping_quality() {
        assert_eq!(MappingQuality::try_from(0), Ok(MappingQuality(0)));
        assert_eq!(MappingQuality::try_from(8), Ok(MappingQuality(8)));
        assert_eq!(MappingQuality::try_from(13), Ok(MappingQuality(13)));
        assert_eq!(MappingQuality::try_from(144), Ok(MappingQuality(144)));
        assert_eq!(MappingQuality::try_from(255), Err(ParseError::Missing));
    }

    #[test]
    fn test_from_mapping_quality_for_u8() {
        assert_eq!(u8::from(MappingQuality(0)), 0);
        assert_eq!(u8::from(MappingQuality(8)), 8);
        assert_eq!(u8::from(MappingQuality(13)), 13);
        assert_eq!(u8::from(MappingQuality(144)), 144);
    }
}
