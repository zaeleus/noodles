use std::{error, fmt};

use noodles_core::{Position, position};

use crate::binning_index::index::header::format::CoordinateSystem;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Parse(position::ParseError),
    /// The position is invalid.
    Invalid(position::TryFromIntError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Parse(e) => Some(e),
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(_) => write!(f, "invalid input"),
            Self::Invalid(_) => write!(f, "invalid input"),
        }
    }
}

pub(super) fn parse_start_position(
    s: &str,
    coordinate_system: CoordinateSystem,
) -> Result<Position, ParseError> {
    match coordinate_system {
        CoordinateSystem::Gff => s.parse().map_err(ParseError::Parse),
        CoordinateSystem::Bed => s
            .parse::<usize>()
            .map_err(ParseError::Parse)
            .and_then(|n| Position::try_from(n + 1).map_err(ParseError::Invalid)),
    }
}
