mod bounds;

use std::io;

use lexical_core::FromLexical;
use noodles_core::Position;

pub(crate) use self::bounds::Bounds;

#[derive(Clone, Eq, PartialEq)]
pub(crate) struct Fields {
    pub(crate) buf: Vec<u8>,
    pub(crate) bounds: Bounds,
}

impl Fields {
    pub(super) fn reference_sequence_name(&self) -> &[u8] {
        &self.buf[self.bounds.reference_sequence_name_range()]
    }

    pub(super) fn feature_start(&self) -> io::Result<Position> {
        let src = &self.buf[self.bounds.feature_start_range()];

        parse_int::<usize>(src).and_then(|n| {
            n.checked_add(1)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "attempt to add with overflow")
                })
                .and_then(|m| {
                    Position::try_from(m).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                })
        })
    }

    pub(super) fn feature_end(&self) -> io::Result<Position> {
        let src = &self.buf[self.bounds.feature_end_range()];

        parse_int::<usize>(src).and_then(|n| {
            Position::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }

    pub(super) fn name(&self) -> Option<&[u8]> {
        const NAME_INDEX: usize = 0;
        self.get(NAME_INDEX)
    }

    pub(super) fn score(&self) -> Option<&[u8]> {
        const SCORE_INDEX: usize = 1;
        self.get(SCORE_INDEX)
    }

    pub(super) fn strand(&self) -> Option<&[u8]> {
        const STRAND_INDEX: usize = 2;
        self.get(STRAND_INDEX)
    }

    pub(super) fn get(&self, i: usize) -> Option<&[u8]> {
        self.bounds.get(i).map(|range| &self.buf[range])
    }
}

impl Default for Fields {
    fn default() -> Self {
        Self {
            buf: Vec::from(*b"sq001"),
            bounds: Bounds::default(),
        }
    }
}

fn parse_int<N: FromLexical>(buf: &[u8]) -> io::Result<N> {
    lexical_core::parse(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}