use std::io;

use noodles_core as core;

/// A raw BAM record position.
#[derive(Debug, Eq, PartialEq)]
pub struct Position(i32);

impl Position {
    pub(super) fn new(n: i32) -> Self {
        Self(n)
    }
}

impl TryFrom<Position> for core::Position {
    type Error = io::Error;

    fn try_from(position: Position) -> Result<Self, Self::Error> {
        usize::try_from(position.0)
            .map(|n| n + 1)
            .and_then(Self::try_from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}
