//! SAM record mapping quality.

use std::ops::Deref;

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
    fn test_from_mapping_quality_for_u8() {
        assert_eq!(u8::from(MappingQuality::from(0)), 0);
        assert_eq!(u8::from(MappingQuality::from(8)), 8);
        assert_eq!(u8::from(MappingQuality::from(13)), 13);
        assert_eq!(u8::from(MappingQuality::from(144)), 144);
        assert_eq!(u8::from(MappingQuality::from(255)), 255);
    }
}
