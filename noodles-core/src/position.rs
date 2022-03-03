use std::num::NonZeroUsize;

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
}

impl From<Position> for usize {
    fn from(position: Position) -> Self {
        position.0.get()
    }
}
