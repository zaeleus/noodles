mod bounds;

use std::io;

use noodles_core::Position;

use self::bounds::Bounds;

#[derive(Clone, Eq, PartialEq)]
pub(super) struct Fields<'l> {
    src: &'l str,
    bounds: Bounds,
}

impl<'l> Fields<'l> {
    pub(super) fn try_new(src: &'l str) -> io::Result<Self> {
        Bounds::index(src).map(|bounds| Self { src, bounds })
    }

    pub fn reference_sequence_name(&self) -> &str {
        &self.src[self.bounds.reference_sequence_name_range()]
    }

    pub fn source(&self) -> &str {
        &self.src[self.bounds.source_range()]
    }

    pub fn ty(&self) -> &str {
        &self.src[self.bounds.type_range()]
    }

    pub fn start(&self) -> io::Result<Position> {
        let src = &self.src[self.bounds.start_range()];
        parse_position(src)
    }

    pub fn end(&self) -> io::Result<Position> {
        let src = &self.src[self.bounds.end_range()];
        parse_position(src)
    }

    pub fn score(&self) -> &str {
        &self.src[self.bounds.score_range()]
    }

    pub fn strand(&self) -> &str {
        &self.src[self.bounds.strand_range()]
    }

    pub fn phase(&self) -> &str {
        &self.src[self.bounds.phase_range()]
    }

    pub fn attributes(&self) -> &str {
        &self.src[self.bounds.attributes_range()]
    }
}

fn parse_position(s: &str) -> io::Result<Position> {
    s.parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
