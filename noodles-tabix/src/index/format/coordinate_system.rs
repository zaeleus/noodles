//! Tabix index format coordinate system.

use std::{convert::TryFrom, error, fmt};

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
    fn from(coordinate_system: CoordinateSystem) -> u16 {
        coordinate_system as u16
    }
}
