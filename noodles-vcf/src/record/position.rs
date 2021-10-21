//! VCF record position.

use std::{error, fmt, num, str::FromStr};

const MIN: i32 = 0;

/// A VCF record position.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Position(i32);

impl From<Position> for i32 {
    fn from(position: Position) -> Self {
        position.0
    }
}

/// An error returned when a raw VCF record position fails to parse.
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

/// An error returned when a raw VCF record position fails to convert.
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
        if n >= MIN {
            Ok(Self(n))
        } else {
            Err(TryFromIntError(n))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_position_for_i32() -> Result<(), TryFromIntError> {
        assert_eq!(Position::try_from(8).map(i32::from)?, 8);
        assert_eq!(Position::try_from(13).map(i32::from)?, 13);
        Ok(())
    }

    #[test]
    fn test_from_str() {
        assert_eq!("13".parse(), Ok(Position(13)));

        assert!(matches!("".parse::<Position>(), Err(ParseError::Parse(_))));
        assert!(matches!(
            "noodles".parse::<Position>(),
            Err(ParseError::Parse(_))
        ));

        assert_eq!(
            "-1".parse::<Position>(),
            Err(ParseError::Invalid(TryFromIntError(-1)))
        );
    }

    #[test]
    fn test_try_from_i32_for_position() {
        assert_eq!(Position::try_from(0), Ok(Position(0)));
        assert_eq!(Position::try_from(13), Ok(Position(13)));

        assert_eq!(Position::try_from(-1), Err(TryFromIntError(-1)));
    }
}
