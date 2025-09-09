//! 1-based position.

mod sequence_index;

pub use self::sequence_index::SequenceIndex;

use std::{
    fmt,
    num::{self, NonZero},
    str::FromStr,
};

/// A 1-based position.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct Position(NonZero<usize>);

impl Position {
    /// The minimum value of a position.
    pub const MIN: Self = Self(NonZero::<usize>::MIN);

    /// The maximum value of a position.
    pub const MAX: Self = Self(NonZero::<usize>::MAX);

    /// Creates a position if the given value is not zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// assert!(Position::new(8).is_some());
    /// assert!(Position::new(0).is_none());
    /// ```
    pub const fn new(n: usize) -> Option<Self> {
        if let Some(m) = NonZero::new(n) {
            Some(Self(m))
        } else {
            None
        }
    }

    /// Returns the inner value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// assert_eq!(Position::MIN.get(), 1);
    /// ```
    pub const fn get(&self) -> usize {
        self.0.get()
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
    pub const fn checked_add(self, other: usize) -> Option<Self> {
        if let Some(n) = self.0.checked_add(other) {
            Some(Self(n))
        } else {
            None
        }
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
        NonZero::try_from(n).map(Position)
    }
}

impl From<Position> for usize {
    fn from(position: Position) -> Self {
        position.0.get()
    }
}
