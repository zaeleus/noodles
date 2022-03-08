//! 1-based position.

mod sequence_index;

pub use self::sequence_index::SequenceIndex;

use std::{
    fmt,
    num::{self, NonZeroUsize},
    str::FromStr,
};

/// A 1-based position.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct Position(NonZeroUsize);

impl Position {
    /// Creates a position if the given value is not zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// assert!(Position::new(8).is_some());
    /// assert!(Position::new(0).is_none());
    /// ```
    pub fn new(n: usize) -> Option<Self> {
        NonZeroUsize::new(n).map(Self)
    }

    /// Adds an unsigned integer to a 1-based position.
    ///
    /// This returns `None` if the operation overflowed.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// let position = Position::try_from(8)?;
    /// assert_eq!(position.checked_add(5), Position::new(13));
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn checked_add(self, other: usize) -> Option<Self> {
        usize::from(self).checked_add(other).and_then(Self::new)
    }
}

impl fmt::Display for Position {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// An error returned when a position fails to parse.
pub type ParseError = num::ParseIntError;

impl FromStr for Position {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.parse().map(Self)
    }
}

/// An error returned when a raw position fails to convert.
pub type TryFromIntError = num::TryFromIntError;

impl TryFrom<usize> for Position {
    type Error = TryFromIntError;

    fn try_from(n: usize) -> Result<Self, Self::Error> {
        NonZeroUsize::try_from(n).map(Position)
    }
}

impl From<Position> for usize {
    fn from(position: Position) -> Self {
        position.0.get()
    }
}
