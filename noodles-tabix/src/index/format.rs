//! Tabix index format and coordinate system.

pub mod coordinate_system;

pub use self::coordinate_system::CoordinateSystem;

use std::{convert::TryFrom, error, fmt};

const COORDINATE_SYSTEM_SHIFT: usize = 16;
const FORMAT_MASK: i32 = 0xffff;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    Generic(CoordinateSystem),
    Sam,
    Vcf,
}

impl Format {
    pub fn coordinate_system(&self) -> CoordinateSystem {
        match self {
            Self::Generic(coordinate_system) => *coordinate_system,
            Self::Sam | Self::Vcf => CoordinateSystem::Gff,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TryFromIntError {
    InvalidCoordinateSystem(coordinate_system::TryFromIntError),
    InvalidKind(u16),
}

impl error::Error for TryFromIntError {}

impl fmt::Display for TryFromIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid format: ")?;

        match self {
            Self::InvalidCoordinateSystem(e) => write!(f, "{}", e),
            Self::InvalidKind(n) => write!(f, "invalid kind: expected 0..=2, got {}", n),
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
            Format::Generic(coordinate_system) => i32::from(u16::from(coordinate_system)) << 16,
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
