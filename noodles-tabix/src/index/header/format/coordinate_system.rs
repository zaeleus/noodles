//! Tabix index format coordinate system.

use std::{error, fmt};

/// A tabix index format coordinate system.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CoordinateSystem {
    /// GFF coordinates: 1-based [start, end]
    Gff,
    /// BED coordinates: 0-based [start, end)
    Bed,
}

/// An error returned when a raw coordinate system fails to convert.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TryFromIntError(u16);

impl error::Error for TryFromIntError {}

impl fmt::Display for TryFromIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid coordinate system: expected 0 or 1, got {}",
            self.0
        )
    }
}

impl TryFrom<u16> for CoordinateSystem {
    type Error = TryFromIntError;

    fn try_from(n: u16) -> Result<Self, Self::Error> {
        match n {
            0 => Ok(Self::Gff),
            1 => Ok(Self::Bed),
            _ => Err(TryFromIntError(n)),
        }
    }
}

impl From<CoordinateSystem> for u16 {
    fn from(coordinate_system: CoordinateSystem) -> Self {
        coordinate_system as u16
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_u16_for_coordinate_system() {
        assert_eq!(CoordinateSystem::try_from(0), Ok(CoordinateSystem::Gff));
        assert_eq!(CoordinateSystem::try_from(1), Ok(CoordinateSystem::Bed));
        assert_eq!(CoordinateSystem::try_from(2), Err(TryFromIntError(2)));
    }

    #[test]
    fn test_from_coordinate_system_for_u16() {
        assert_eq!(u16::from(CoordinateSystem::Gff), 0);
        assert_eq!(u16::from(CoordinateSystem::Bed), 1);
    }
}
