//! Genomic region interval.

use std::{
    error, fmt,
    ops::{RangeBounds, RangeFrom, RangeFull, RangeInclusive, RangeToInclusive},
    str::FromStr,
};

use crate::{position, Position};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Bound {
    Included(Position),
    Unbounded,
}

/// A closed interval.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Interval {
    start: Bound,
    end: Bound,
}

impl Interval {
    /// Creates a closed interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::ops::{Bound, RangeBounds};
    /// use noodles_core::{region::Interval, Position};
    ///
    /// let start = Position::try_from(8)?;
    /// let end = Position::try_from(13)?;
    ///
    /// let interval = Interval::new(start, end);
    /// assert_eq!(interval.start_bound().cloned(), Bound::Included(start));
    /// assert_eq!(interval.end_bound().cloned(), Bound::Included(end));
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn new(start: Position, end: Position) -> Self {
        Self {
            start: Bound::Included(start),
            end: Bound::Included(end),
        }
    }

    /// Resolves the start position and returns it as an inclusive start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::{region::Interval, Position};
    /// let start = Position::try_from(8)?;
    /// let end = Position::try_from(13)?;
    /// let interval = Interval::new(start, end);
    /// assert_eq!(interval.start(), start);
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn start(&self) -> Position {
        match self.start {
            Bound::Included(start) => start,
            Bound::Unbounded => Position::MIN,
        }
    }

    /// Resolves the end position and returns it as an inclusive end.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::{region::Interval, Position};
    /// let start = Position::try_from(8)?;
    /// let end = Position::try_from(13)?;
    /// let interval = Interval::new(start, end);
    /// assert_eq!(interval.start(), start);
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn end(&self) -> Position {
        match self.end {
            Bound::Included(end) => end,
            Bound::Unbounded => Position::MAX,
        }
    }
}

impl fmt::Display for Interval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (self.start, self.end) {
            (Bound::Unbounded, Bound::Unbounded) => Ok(()),
            (Bound::Unbounded, Bound::Included(e)) => write!(f, "{}-{}", Position::MIN, e),
            (Bound::Included(s), Bound::Unbounded) => s.fmt(f),
            (Bound::Included(s), Bound::Included(e)) => write!(f, "{}-{}", s, e),
        }
    }
}

impl RangeBounds<Position> for Interval {
    fn start_bound(&self) -> std::ops::Bound<&Position> {
        bound_to_std_ops_bound(&self.start)
    }

    fn end_bound(&self) -> std::ops::Bound<&Position> {
        bound_to_std_ops_bound(&self.end)
    }
}

fn bound_to_std_ops_bound(bound: &Bound) -> std::ops::Bound<&Position> {
    match bound {
        Bound::Included(ref value) => std::ops::Bound::Included(value),
        Bound::Unbounded => std::ops::Bound::Unbounded,
    }
}

/// An error returned when an interval fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The start position is invalid.
    InvalidStartPosition(position::ParseError),
    /// The end position is invalid.
    InvalidEndPosition(position::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidStartPosition(e) => write!(f, "invalid start position: {}", e),
            Self::InvalidEndPosition(e) => write!(f, "invalid end position: {}", e),
        }
    }
}

impl FromStr for Interval {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Ok(Self::from(..));
        }

        let mut components = s.splitn(2, '-');

        let start = match components.next() {
            Some(t) => t
                .parse()
                .map(Bound::Included)
                .map_err(ParseError::InvalidStartPosition)?,
            None => Bound::Unbounded,
        };

        let end = match components.next() {
            Some(t) => t
                .parse()
                .map(Bound::Included)
                .map_err(ParseError::InvalidEndPosition)?,
            None => Bound::Unbounded,
        };

        Ok(Self { start, end })
    }
}

impl From<RangeFrom<Position>> for Interval {
    fn from(range: RangeFrom<Position>) -> Self {
        Self {
            start: Bound::Included(range.start),
            end: Bound::Unbounded,
        }
    }
}

impl From<RangeFull> for Interval {
    fn from(_: RangeFull) -> Self {
        Self {
            start: Bound::Unbounded,
            end: Bound::Unbounded,
        }
    }
}

impl From<RangeInclusive<Position>> for Interval {
    fn from(range: RangeInclusive<Position>) -> Self {
        let (start, end) = range.into_inner();

        Self {
            start: Bound::Included(start),
            end: Bound::Included(end),
        }
    }
}

impl From<RangeToInclusive<Position>> for Interval {
    fn from(range: RangeToInclusive<Position>) -> Self {
        Self {
            start: Bound::Unbounded,
            end: Bound::Included(range.end),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() -> Result<(), crate::position::TryFromIntError> {
        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;

        let interval = Interval {
            start: Bound::Unbounded,
            end: Bound::Unbounded,
        };
        assert_eq!(interval.to_string(), "");

        let interval = Interval {
            start: Bound::Unbounded,
            end: Bound::Included(end),
        };
        assert_eq!(interval.to_string(), "1-13");

        let interval = Interval {
            start: Bound::Included(start),
            end: Bound::Unbounded,
        };
        assert_eq!(interval.to_string(), "8");

        assert_eq!(Interval::new(start, end).to_string(), "8-13");

        Ok(())
    }

    #[test]
    fn test_parse() -> Result<(), crate::position::TryFromIntError> {
        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;

        assert_eq!(
            "".parse(),
            Ok(Interval {
                start: Bound::Unbounded,
                end: Bound::Unbounded
            })
        );

        assert_eq!(
            "8".parse(),
            Ok(Interval {
                start: Bound::Included(start),
                end: Bound::Unbounded
            })
        );

        assert_eq!(
            "8-13".parse(),
            Ok(Interval {
                start: Bound::Included(start),
                end: Bound::Included(end),
            })
        );

        assert!(matches!(
            "x".parse::<Interval>(),
            Err(ParseError::InvalidStartPosition(_))
        ));

        assert!(matches!(
            "1-x".parse::<Interval>(),
            Err(ParseError::InvalidEndPosition(_))
        ));

        Ok(())
    }
}
