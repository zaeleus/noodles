//! Tabix index header format and coordinate system.

pub mod coordinate_system;

pub use self::coordinate_system::CoordinateSystem;

use std::{error, fmt};

const COORDINATE_SYSTEM_SHIFT: usize = 16;
const FORMAT_MASK: i32 = 0xffff;

/// A tabix index format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    /// A generic format with a defined coordinate system.
    Generic(CoordinateSystem),
    /// The SAM (Sequence Alignment/Map) format.
    Sam,
    /// The VCF (Variant Call Format) format.
    Vcf,
}

impl Format {
    /// Returns the coordinate system of the format.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::index::header::{format::CoordinateSystem, Format};
    ///
    /// let format = Format::Generic(CoordinateSystem::Bed);
    /// assert_eq!(format.coordinate_system(), CoordinateSystem::Bed);
    ///
    /// assert_eq!(Format::Sam.coordinate_system(), CoordinateSystem::Gff);
    /// assert_eq!(Format::Vcf.coordinate_system(), CoordinateSystem::Gff);
    /// ```
    pub fn coordinate_system(&self) -> CoordinateSystem {
        match self {
            Self::Generic(coordinate_system) => *coordinate_system,
            Self::Sam | Self::Vcf => CoordinateSystem::Gff,
        }
    }
}

/// An error returned when a raw format fails to convert.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TryFromIntError {
    /// The coordinate system is invalid.
    InvalidCoordinateSystem(coordinate_system::TryFromIntError),
    /// The kind is invalid.
    InvalidKind(u16),
}

impl error::Error for TryFromIntError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidCoordinateSystem(e) => Some(e),
            Self::InvalidKind(_) => None,
        }
    }
}

impl fmt::Display for TryFromIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidCoordinateSystem(_) => f.write_str("invalid coordinate system"),
            Self::InvalidKind(n) => write!(f, "invalid kind: expected 0..=2, got {n}"),
        }
    }
}

impl TryFrom<i32> for Format {
    type Error = TryFromIntError;

    fn try_from(n: i32) -> Result<Self, Self::Error> {
        let raw_kind = (n & FORMAT_MASK) as u16;

        match raw_kind {
            0 => {
                let raw_format = (n >> COORDINATE_SYSTEM_SHIFT) as u16;
                CoordinateSystem::try_from(raw_format)
                    .map(Self::Generic)
                    .map_err(TryFromIntError::InvalidCoordinateSystem)
            }
            1 => Ok(Self::Sam),
            2 => Ok(Self::Vcf),
            _ => Err(TryFromIntError::InvalidKind(raw_kind)),
        }
    }
}

impl From<Format> for i32 {
    fn from(format: Format) -> Self {
        match format {
            Format::Generic(coordinate_system) => {
                i32::from(u16::from(coordinate_system)) << COORDINATE_SYSTEM_SHIFT
            }
            Format::Sam => 1,
            Format::Vcf => 2,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_i32_for_format() {
        assert_eq!(
            Format::try_from(0x010000),
            Ok(Format::Generic(CoordinateSystem::Bed))
        );

        assert_eq!(
            Format::try_from(0),
            Ok(Format::Generic(CoordinateSystem::Gff))
        );

        assert_eq!(Format::try_from(1), Ok(Format::Sam));
        assert_eq!(Format::try_from(2), Ok(Format::Vcf));

        assert_eq!(Format::try_from(3), Err(TryFromIntError::InvalidKind(3)));
    }

    #[test]
    fn test_from_format_for_i32() {
        assert_eq!(i32::from(Format::Generic(CoordinateSystem::Bed)), 0x010000);
        assert_eq!(i32::from(Format::Generic(CoordinateSystem::Gff)), 0);
        assert_eq!(i32::from(Format::Sam), 1);
        assert_eq!(i32::from(Format::Vcf), 2);
    }
}
