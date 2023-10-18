use std::io;

use noodles_core as core;

/// Raw SAM record position.
#[derive(Debug, Eq, PartialEq)]
pub struct Position<'a>(&'a [u8]);

impl<'a> Position<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
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
        lexical_core::parse(raw_position.as_ref())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|n: usize| {
                core::Position::try_from(n)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }
}
