//! SAM record mapping quality.

use std::{error, fmt, num, ops::Deref, str::FromStr};

/// The raw value of a missing mapping quality.
pub const MISSING: u8 = 255;

/// A SAM record mapping quality.
///
/// Mapping quality ranges from 0 to 254 (inclusive), where higher is better.
///
/// The value 255 is reserved as a marker for a missing mapping quality.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct MappingQuality(Option<u8>);

impl From<u8> for MappingQuality {
    fn from(n: u8) -> Self {
        if n == MISSING {
            Self(None)
        } else {
            Self(Some(n))
        }
    }
}

/// An error returned when a raw SAM record mapping quality fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input failed to parse as an integer.
    Parse(num::ParseIntError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(e) => write!(f, "parse error: {}", e),
        }
    }
}

impl FromStr for MappingQuality {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let n: u8 = s.parse().map_err(ParseError::Parse)?;
        Ok(Self::from(n))
    }
}

impl From<MappingQuality> for u8 {
    fn from(mapping_quality: MappingQuality) -> Self {
        mapping_quality.unwrap_or(MISSING)
    }
}

impl Deref for MappingQuality {
    type Target = Option<u8>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_u8_for_mapping_quality() {
        assert_eq!(*MappingQuality::from(0), Some(0));
        assert_eq!(*MappingQuality::from(8), Some(8));
        assert_eq!(*MappingQuality::from(13), Some(13));
        assert_eq!(*MappingQuality::from(144), Some(144));
        assert_eq!(*MappingQuality::from(255), None);
    }

    #[test]
    fn test_from_str() {
        assert_eq!("0".parse(), Ok(MappingQuality::from(0)));
        assert_eq!("8".parse(), Ok(MappingQuality::from(8)));
        assert_eq!("13".parse(), Ok(MappingQuality::from(13)));
        assert_eq!("144".parse(), Ok(MappingQuality::from(144)));
        assert_eq!("255".parse(), Ok(MappingQuality::from(255)));

        assert!(matches!(
            "".parse::<MappingQuality>(),
            Err(ParseError::Parse(_))
        ));
        assert!(matches!(
            "256".parse::<MappingQuality>(),
            Err(ParseError::Parse(_))
        ));
    }

    #[test]
    fn test_from_mapping_quality_for_u8() {
        assert_eq!(u8::from(MappingQuality::from(0)), 0);
        assert_eq!(u8::from(MappingQuality::from(8)), 8);
        assert_eq!(u8::from(MappingQuality::from(13)), 13);
        assert_eq!(u8::from(MappingQuality::from(144)), 144);
        assert_eq!(u8::from(MappingQuality::from(255)), 255);
    }
}
