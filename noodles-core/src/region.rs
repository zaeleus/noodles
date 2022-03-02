//! Genomic region.

mod mapped;

pub use self::mapped::{Interval, Mapped};

use std::{
    error, fmt, num,
    ops::{Bound, RangeBounds},
    str::FromStr,
};

static UNMAPPED_NAME: &str = "*";
static ALL_NAME: &str = ".";

/// A genomic region.
///
/// Genomic regions can either be mapped to a reference sequence, unmapped (*), or an inclusion of
/// all records (.).
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Region {
    /// A mapped region.
    Mapped(Mapped),
    /// An unmapped region.
    Unmapped,
    /// All records.
    All,
}

impl Region {
    /// Creates a new mapped region.
    ///
    /// Positions are assumed to be 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::ops::Bound;
    /// use noodles_core::Region;
    /// let region = Region::mapped("sq0", 1..=5);
    /// assert!(matches!(region, Region::Mapped(_)));
    /// ```
    pub fn mapped<I, B>(name: I, interval: B) -> Self
    where
        I: Into<String>,
        B: RangeBounds<i32>,
    {
        Self::Mapped(Mapped::new(name, interval))
    }

    /// Returns the reference name of the region.
    ///
    /// If the region is unmapped, this returns "*". If the region represents
    /// all records, this returns ".".
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Region;
    ///
    /// let region = Region::mapped("sq0", 1..=5);
    /// assert_eq!(region.name(), "sq0");
    ///
    /// assert_eq!(Region::Unmapped.name(), "*");
    /// assert_eq!(Region::All.name(), ".");
    /// ```
    pub fn name(&self) -> &str {
        match self {
            Self::Mapped(m) => m.name(),
            Self::Unmapped => UNMAPPED_NAME,
            Self::All => ALL_NAME,
        }
    }

    /// Returns whether the region is mapped.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Region;
    /// assert!(Region::mapped("sq0", 5..=8).is_mapped());
    /// assert!(!Region::Unmapped.is_mapped());
    /// ```
    pub fn is_mapped(&self) -> bool {
        matches!(self, Self::Mapped(_))
    }

    /// Returns the region as a mapped region if it is mapped.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Region;
    /// assert!(Region::mapped("sq0", 5..=8).as_mapped().is_some());
    /// assert!(Region::Unmapped.as_mapped().is_none());
    /// ```
    pub fn as_mapped(&self) -> Option<&Mapped> {
        match self {
            Self::Mapped(m) => Some(m),
            _ => None,
        }
    }
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Mapped(m) => write!(f, "{}", m),
            Self::Unmapped => write!(f, "{}", UNMAPPED_NAME),
            Self::All => write!(f, "{}", ALL_NAME),
        }
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
    /// The start position is invalid.
    InvalidStartPosition(num::ParseIntError),
    /// The end position is invalid.
    InvalidEndPosition(num::ParseIntError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Ambiguous => f.write_str("ambiguous input"),
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidStartPosition(e) => write!(f, "invalid start position: {}", e),
            Self::InvalidEndPosition(e) => write!(f, "invalid end position: {}", e),
        }
    }
}

impl FromStr for Region {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        } else if s == UNMAPPED_NAME {
            return Ok(Self::Unmapped);
        } else if s == ALL_NAME {
            return Ok(Self::All);
        }

        if let Some((name, suffix)) = s.rsplit_once(':') {
            if let Ok(interval) = parse_interval(suffix) {
                Ok(Self::mapped(name, interval))
            } else {
                Err(ParseError::Invalid)
            }
        } else {
            Ok(Self::mapped(s, ..))
        }
    }
}

fn parse_interval(s: &str) -> Result<Interval, ParseError> {
    if s.is_empty() {
        return Ok((Bound::Unbounded, Bound::Unbounded));
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

    Ok((start, end))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Region::mapped("sq0", ..).to_string(), "sq0");
        assert_eq!(Region::mapped("sq0", 3..).to_string(), "sq0:3");
        assert_eq!(Region::mapped("sq0", 3..=5).to_string(), "sq0:3-5");
        assert_eq!(Region::Unmapped.to_string(), "*");
        assert_eq!(Region::All.to_string(), ".");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("*".parse(), Ok(Region::Unmapped));
        assert_eq!(".".parse(), Ok(Region::All));

        assert_eq!("sq0".parse(), Ok(Region::mapped("sq0", ..)));
        assert_eq!("sq1:".parse(), Ok(Region::mapped("sq1", ..)));
        assert_eq!("sq2:5".parse(), Ok(Region::mapped("sq2", 5..)));
        assert_eq!("sq3:5-8".parse(), Ok(Region::mapped("sq3", 5..=8)));

        assert_eq!("".parse::<Region>(), Err(ParseError::Empty));
    }
}
