use std::io;

use noodles_core as core;

/// An alignment record 1-based position.
pub trait Position {
    /// Converts a position to a `usize`.
    fn try_to_usize(&self) -> io::Result<usize>;
}

impl TryFrom<&dyn Position> for usize {
    type Error = io::Error;

    fn try_from(raw_position: &dyn Position) -> Result<Self, Self::Error> {
        raw_position.try_to_usize()
    }
}

impl TryFrom<&dyn Position> for core::Position {
    type Error = io::Error;

    fn try_from(raw_position: &dyn Position) -> Result<Self, Self::Error> {
        usize::try_from(raw_position).and_then(|n| {
            Self::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }
}
