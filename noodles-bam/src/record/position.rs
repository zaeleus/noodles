use std::io;

use noodles_core as core;
use noodles_sam::{self as sam, alignment::record::fields::Position as _};

/// A raw BAM record position.
#[derive(Debug, Eq, PartialEq)]
pub struct Position(i32);

impl Position {
    pub(super) fn new(n: i32) -> Self {
        Self(n)
    }
}

impl sam::alignment::record::fields::Position for Position {
    fn try_to_usize(&self) -> io::Result<usize> {
        usize::try_from(self.0)
            .map(|n| n + 1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}

impl TryFrom<Position> for core::Position {
    type Error = io::Error;

    fn try_from(position: Position) -> Result<Self, Self::Error> {
        position.try_to_usize().and_then(|n| {
            Self::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }
}
