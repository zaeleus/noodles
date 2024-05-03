//! Genomic region.

pub mod interval;

use bstr::{BStr, BString};

pub use self::interval::Interval;

use std::{
    error, fmt,
    ops::{Bound, RangeBounds},
    str::FromStr,
};

use super::Position;

/// A genomic region.
///
/// A genomic region describes a region on a reference sequence. It consists of reference sequence
/// name and an interval.
///
/// They are represented in text as `reference-sequence-name[:start[-end]]`, where the start and
/// end positions are 1-based, inclusive. If no end position is given, it is assumed to span from
/// the start to the end of the reference sequence. If no interval is given, it is assumed to span
/// the entirety of the reference sequence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Region {
    name: BString,
    interval: Interval,
}

impl Region {
    /// Creates a region.
    ///
    /// Positions are assumed to be 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::{Region, Position};
    ///
    /// let start = Position::try_from(5)?;
    /// let end = Position::try_from(8)?;
    /// let region = Region::new("sq0", start..=end);
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn new<N, I>(name: N, interval: I) -> Self
    where
        N: Into<BString>,
        I: Into<Interval>,
    {
        Self {
            name: name.into(),
            interval: interval.into(),
        }
    }

    /// Returns the reference name of the region.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::{Position, Region};
    ///
    /// let start = Position::try_from(5)?;
    /// let end = Position::try_from(8)?;
    /// let region = Region::new("sq0", start..=end);
    ///
    /// assert_eq!(region.name(), &b"sq0"[..]);
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn name(&self) -> &BStr {
        self.name.as_ref()
    }

    /// Returns the start position of the region (1-based).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::ops::Bound;
    /// use noodles_core::{Position, Region};
    ///
    /// let start = Position::try_from(5)?;
    /// let region = Region::new("sq0", start..);
    ///
    /// assert_eq!(region.start(), Bound::Included(start));
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn start(&self) -> Bound<Position> {
        self.interval.start_bound().cloned()
    }

    /// Returns the end position of the region (1-based).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::ops::Bound;
    /// use noodles_core::{Position, Region};
    ///
    /// let end = Position::try_from(8)?;
    /// let region = Region::new("sq0", ..=end);
    ///
    /// assert_eq!(region.end(), Bound::Included(end));
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn end(&self) -> Bound<Position> {
        self.interval.end_bound().cloned()
    }

    /// Returns the start and end positions as an interval.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::ops::Bound;
    /// use noodles_core::{region::Interval, Position, Region};
    ///
    /// let start = Position::try_from(5)?;
    /// let end = Position::try_from(8)?;
    /// let region = Region::new("sq0", start..=end);
    ///
    /// assert_eq!(region.interval(), Interval::from(start..=end));
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn interval(&self) -> Interval {
        self.interval
    }
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name())?;

        match (self.interval.start_bound(), self.interval.end_bound()) {
            (Bound::Unbounded, Bound::Unbounded) => {}
            (_, _) => write!(f, ":{}", self.interval)?,
        }

        Ok(())
    }
}

/// An error returned when a genomic region fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is ambiguous.
    Ambiguous,
    /// The input is invalid.
    Invalid,
    /// The interval is invalid.
    InvalidInterval(interval::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidInterval(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Ambiguous => f.write_str("ambiguous input"),
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidInterval(_) => f.write_str("invalid interval"),
        }
    }
}

impl FromStr for Region {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        if let Some((name, suffix)) = s.rsplit_once(':') {
            let interval: Interval = suffix.parse().map_err(ParseError::InvalidInterval)?;
            Ok(Self::new(name, interval))
        } else {
            Ok(Self::new(s, ..))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() -> Result<(), crate::position::TryFromIntError> {
        let start = Position::try_from(5)?;
        let end = Position::try_from(8)?;

        assert_eq!(Region::new("sq0", ..).to_string(), "sq0");
        assert_eq!(Region::new("sq0", ..=end).to_string(), "sq0:1-8");
        assert_eq!(Region::new("sq0", start..).to_string(), "sq0:5");
        assert_eq!(Region::new("sq0", start..=end).to_string(), "sq0:5-8");

        Ok(())
    }

    #[test]
    fn test_from_str() -> Result<(), crate::position::TryFromIntError> {
        assert_eq!("sq0".parse(), Ok(Region::new("sq0", ..)));
        assert_eq!("sq1:".parse(), Ok(Region::new("sq1", ..)));

        let start = Position::try_from(5)?;
        assert_eq!("sq2:5".parse(), Ok(Region::new("sq2", start..)));

        let end = Position::try_from(8)?;
        assert_eq!("sq3:5-8".parse(), Ok(Region::new("sq3", start..=end)));

        assert_eq!("".parse::<Region>(), Err(ParseError::Empty));

        Ok(())
    }
}
