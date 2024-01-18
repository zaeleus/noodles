use std::io;

use noodles_core as core;

use crate::alignment::record::fields::Position as _;

/// Raw SAM record position.
#[derive(Debug, Eq, PartialEq)]
pub struct Position<'a>(&'a [u8]);

impl<'a> Position<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }
}

impl<'a> crate::alignment::record::fields::Position for Position<'a> {
    fn try_to_usize(&self) -> io::Result<usize> {
        lexical_core::parse(self.as_ref())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}

impl<'a> AsRef<[u8]> for Position<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryFrom<Position<'a>> for core::Position {
    type Error = io::Error;

    fn try_from(raw_position: Position<'a>) -> Result<Self, Self::Error> {
        raw_position.try_to_usize().and_then(|n: usize| {
            core::Position::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }
}
