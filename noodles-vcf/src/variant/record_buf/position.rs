//! VCF record position.

use std::{cmp::Ordering, fmt, num, str::FromStr};

use noodles_core as core;

/// A VCF record position.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Position(usize);

impl fmt::Display for Position {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// An error returned when a raw VCF record position fails to parse.
pub type ParseError = num::ParseIntError;

impl FromStr for Position {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.parse().map(Self)
    }
}

impl From<usize> for Position {
    fn from(n: usize) -> Self {
        Self(n)
    }
}

impl From<core::Position> for Position {
    fn from(position: core::Position) -> Self {
        Self::from(usize::from(position))
    }
}

impl From<Position> for usize {
    fn from(position: Position) -> Self {
        position.0
    }
}

impl PartialEq<core::Position> for Position {
    fn eq(&self, other: &core::Position) -> bool {
        self.0.eq(&usize::from(*other))
    }
}

impl PartialOrd<core::Position> for Position {
    fn partial_cmp(&self, other: &core::Position) -> Option<Ordering> {
        if self.0 == 0 {
            None
        } else {
            self.0.partial_cmp(&usize::from(*other))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Position(0).to_string(), "0");
        assert_eq!(Position(8).to_string(), "8");
        assert_eq!(Position(13).to_string(), "13");
    }

    #[test]
    fn test_from_str() {
        use std::num::IntErrorKind;

        assert_eq!("0".parse(), Ok(Position(0)));
        assert_eq!("8".parse(), Ok(Position(8)));
        assert_eq!("13".parse(), Ok(Position(13)));

        assert!(matches!("".parse::<Position>(), Err(e) if e.kind() == &IntErrorKind::Empty));
        assert!(
            matches!("ndls".parse::<Position>(), Err(e) if e.kind() == &IntErrorKind::InvalidDigit)
        );
        assert!(
            matches!("-1".parse::<Position>(), Err(e) if e.kind() == &IntErrorKind::InvalidDigit)
        );
    }

    #[test]
    fn test_from_usize_for_position() {
        assert_eq!(Position::from(0), Position(0));
        assert_eq!(Position::from(8), Position(8));
        assert_eq!(Position::from(13), Position(13));
    }

    #[test]
    fn test_from_position_for_usize() {
        assert_eq!(usize::from(Position::from(0)), 0);
        assert_eq!(usize::from(Position::from(8)), 8);
        assert_eq!(usize::from(Position::from(13)), 13);
    }

    #[test]
    fn test_partial_eq_core_position_for_position() {
        let q = core::Position::MIN;

        let p = Position::from(1);
        assert_eq!(p, q);

        let p = Position::from(0);
        assert_ne!(p, q);
    }

    #[allow(clippy::nonminimal_bool)]
    #[test]
    fn test_partial_ord_core_position_for_position() -> Result<(), core::position::TryFromIntError>
    {
        let q = core::Position::try_from(8)?;

        let p = Position::from(7);
        assert!(p < q);
        assert!(p <= q);
        assert!(!(p >= q));
        assert!(!(p > q));

        let p = Position::from(8);
        assert!(!(p < q));
        assert!(p <= q);
        assert!(p >= q);
        assert!(!(p > q));

        let p = Position::from(9);
        assert!(!(p < q));
        assert!(!(p <= q));
        assert!(p >= q);
        assert!(p > q);

        let p = Position::from(0);
        assert!(!(p < q));
        assert!(!(p <= q));
        assert!(!(p >= q));
        assert!(!(p > q));

        Ok(())
    }
}
