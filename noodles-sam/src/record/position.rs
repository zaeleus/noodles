//! SAM record position.

use std::{
    error, fmt,
    num::{self, NonZeroI32},
    ops::Deref,
    str::FromStr,
};

const UNMAPPED: i32 = 0;

/// A SAM record position.
///
/// This represents a 1-based start position on the reference sequence.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Position(Option<NonZeroI32>);

impl From<i32> for Position {
    fn from(n: i32) -> Self {
        if n < UNMAPPED {
            Self(None)
        } else {
            Self(NonZeroI32::new(n))
        }
    }
}

impl From<Position> for i32 {
    fn from(position: Position) -> Self {
        position.map(i32::from).unwrap_or(UNMAPPED)
    }
}

/// An error returned when a raw SAM record position fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input failed to parse as an integer.
    Parse(num::ParseIntError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(e) => write!(f, "{}", e),
        }
    }
}

impl FromStr for Position {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let n: i32 = s.parse().map_err(ParseError::Parse)?;
        Ok(Self::from(n))
    }
}

impl Deref for Position {
    type Target = Option<NonZeroI32>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_i32_for_position() {
        assert_eq!(*Position::from(0), None);
        assert_eq!(*Position::from(13), NonZeroI32::new(13));
    }

    #[test]
    fn test_from_position_for_i32() {
        assert_eq!(i32::from(Position::from(0)), 0);
        assert_eq!(i32::from(Position::from(13)), 13);
    }

    #[test]
    fn test_from_str() {
        assert_eq!("0".parse(), Ok(Position(None)));
        assert_eq!("13".parse(), Ok(Position(NonZeroI32::new(13))));

        assert!(matches!("".parse::<Position>(), Err(ParseError::Parse(_))));
        assert!(matches!(
            "noodles".parse::<Position>(),
            Err(ParseError::Parse(_))
        ));
    }
}
