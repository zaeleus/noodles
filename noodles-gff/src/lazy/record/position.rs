use noodles_core as core;

/// Raw GFF record position.
#[derive(Debug)]
pub struct Position<'a>(&'a str);

impl<'a> Position<'a> {
    pub(super) fn new(buf: &'a str) -> Self {
        Self(buf)
    }
}

impl<'a> AsRef<str> for Position<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

impl<'a> TryFrom<Position<'a>> for core::Position {
    type Error = core::position::ParseError;

    fn try_from(raw_position: Position<'a>) -> Result<Self, Self::Error> {
        raw_position.as_ref().parse()
    }
}
