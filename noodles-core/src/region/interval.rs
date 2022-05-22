//! Genomic region interval.

use std::{
    error, fmt,
    ops::{
        Bound, Range, RangeBounds, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive,
    },
    str::FromStr,
};

use crate::{position, Position};

/// An interval.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Interval {
    start: Bound<Position>,
    end: Bound<Position>,
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
}

impl fmt::Display for Interval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (self.start, self.end) {
            (Bound::Unbounded, Bound::Unbounded) => Ok(()),
            (Bound::Unbounded, Bound::Included(e)) => write!(f, "{}-{}", Position::MIN, e),
            (Bound::Included(s), Bound::Unbounded) => s.fmt(f),
            (Bound::Included(s), Bound::Included(e)) => write!(f, "{}-{}", s, e),
            _ => todo!(),
        }
    }
}

impl RangeBounds<Position> for Interval {
    fn start_bound(&self) -> Bound<&Position> {
        bound_as_ref(&self.start)
    }

    fn end_bound(&self) -> Bound<&Position> {
        bound_as_ref(&self.end)
    }
}

// https://doc.rust-lang.org/nightly/unstable-book/library-features/bound-as-ref.html
fn bound_as_ref<T>(bound: &Bound<T>) -> Bound<&T> {
    match bound {
        Bound::Included(ref value) => Bound::Included(value),
        Bound::Excluded(ref value) => Bound::Excluded(value),
        Bound::Unbounded => Bound::Unbounded,
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

        Ok(Self::from((start, end)))
    }
}

impl From<Range<Position>> for Interval {
    fn from(range: Range<Position>) -> Self {
        Self::from((Bound::Included(range.start), Bound::Excluded(range.end)))
    }
}

impl From<RangeFrom<Position>> for Interval {
    fn from(range: RangeFrom<Position>) -> Self {
        Self::from((Bound::Included(range.start), Bound::Unbounded))
    }
}

impl From<RangeFull> for Interval {
    fn from(_: RangeFull) -> Self {
        Self::from((Bound::Unbounded, Bound::Unbounded))
    }
}

impl From<RangeInclusive<Position>> for Interval {
    fn from(range: RangeInclusive<Position>) -> Self {
        Self::from((range.start_bound().cloned(), range.end_bound().cloned()))
    }
}

impl From<RangeTo<Position>> for Interval {
    fn from(range: RangeTo<Position>) -> Self {
        Self::from((Bound::Unbounded, Bound::Excluded(range.end)))
    }
}

impl From<RangeToInclusive<Position>> for Interval {
    fn from(range: RangeToInclusive<Position>) -> Self {
        Self::from((Bound::Unbounded, Bound::Included(range.end)))
    }
}

impl From<(Bound<Position>, Bound<Position>)> for Interval {
    fn from(bounds: (Bound<Position>, Bound<Position>)) -> Self {
        Self {
            start: bounds.0,
            end: bounds.1,
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
