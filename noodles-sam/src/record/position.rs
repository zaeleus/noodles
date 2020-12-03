//! SAM record position.

use std::{
    convert::TryFrom,
    error, fmt,
    num::{self, NonZeroI32},
    str::FromStr,
};

pub(crate) const UNMAPPED: i32 = 0;

/// A SAM record position.
///
/// This represents a 1-based start position on the reference sequence. The value is guaranteed to
/// be a positive, non-zero integer.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Position(NonZeroI32);

impl From<Position> for i32 {
    fn from(position: Position) -> Self {
        Self::from(position.0)
    }
}

/// An error returned when a raw SAM record position fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input failed to parse as an integer.
    Parse(num::ParseIntError),
    /// The input is invalid.
    Invalid(TryFromIntError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(e) => write!(f, "{}", e),
            Self::Invalid(e) => write!(f, "{}", e),
        }
    }
}

impl FromStr for Position {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let n: i32 = s.parse().map_err(ParseError::Parse)?;
        Self::try_from(n).map_err(ParseError::Invalid)
    }
}

/// An error returned when a raw SAM record position fails to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromIntError(i32);

impl error::Error for TryFromIntError {}

impl fmt::Display for TryFromIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid value: {}", self.0)
    }
}

impl TryFrom<i32> for Position {
    type Error = TryFromIntError;

    fn try_from(n: i32) -> Result<Self, Self::Error> {
        if n < 0 {
            Err(TryFromIntError(n))
        } else {
            NonZeroI32::new(n).map(Self).ok_or(TryFromIntError(n))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_position_for_i32() -> Result<(), TryFromIntError> {
        assert_eq!(i32::from(Position::try_from(8)?), 8);
        assert_eq!(i32::from(Position::try_from(13)?), 13);
        Ok(())
    }

    #[test]
    fn test_from_str() -> Result<(), TryFromIntError> {
        assert_eq!("13".parse(), Ok(Position::try_from(13)?));

        assert!(matches!("".parse::<Position>(), Err(ParseError::Parse(_))));
        assert!(matches!(
            "noodles".parse::<Position>(),
            Err(ParseError::Parse(_))
        ));

        assert_eq!(
            "0".parse::<Position>(),
            Err(ParseError::Invalid(TryFromIntError(0)))
        );
        assert_eq!(
            "-8".parse::<Position>(),
            Err(ParseError::Invalid(TryFromIntError(-8)))
        );

        Ok(())
    }

    #[test]
    fn test_try_from_i32_for_position() -> Result<(), num::TryFromIntError> {
        assert_eq!(
            Position::try_from(13),
            Ok(Position(NonZeroI32::try_from(13)?))
        );

        assert_eq!(Position::try_from(0), Err(TryFromIntError(0)));
        assert_eq!(Position::try_from(-8), Err(TryFromIntError(-8)));

        Ok(())
    }
}
