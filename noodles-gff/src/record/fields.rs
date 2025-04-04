mod bounds;

use std::io;

use bstr::{BStr, ByteSlice};
use noodles_core::Position;

use self::bounds::Bounds;
use super::Attributes;

#[derive(Clone, Eq, PartialEq)]
pub(super) struct Fields<'l> {
    src: &'l [u8],
    bounds: Bounds,
}

impl<'l> Fields<'l> {
    pub(super) fn try_new(src: &'l [u8]) -> io::Result<Self> {
        Bounds::index(src).map(|bounds| Self { src, bounds })
    }

    pub fn reference_sequence_name(&self) -> &BStr {
        self.src[self.bounds.reference_sequence_name_range()].as_bstr()
    }

    pub fn source(&self) -> &BStr {
        self.src[self.bounds.source_range()].as_bstr()
    }

    pub fn ty(&self) -> &BStr {
        self.src[self.bounds.type_range()].as_bstr()
    }

    pub fn start(&self) -> io::Result<Position> {
        let src = &self.src[self.bounds.start_range()];
        parse_position(src)
    }

    pub fn end(&self) -> io::Result<Position> {
        let src = &self.src[self.bounds.end_range()];
        parse_position(src)
    }

    pub fn score(&self) -> &BStr {
        self.src[self.bounds.score_range()].as_bstr()
    }

    pub fn strand(&self) -> &BStr {
        self.src[self.bounds.strand_range()].as_bstr()
    }

    pub fn phase(&self) -> &BStr {
        self.src[self.bounds.phase_range()].as_bstr()
    }

    pub fn attributes(&self) -> Attributes<'_> {
        const MISSING: &[u8] = b".";

        match &self.src[self.bounds.attributes_range()] {
            MISSING => Attributes::new(b""),
            buf => Attributes::new(buf),
        }
    }
}

fn parse_position(src: &[u8]) -> io::Result<Position> {
    lexical_core::parse::<usize>(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|n| {
            Position::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}
