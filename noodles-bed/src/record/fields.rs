mod bounds;

use std::io;

use lexical_core::FromLexical;
use noodles_core::Position;

pub(crate) use self::bounds::Bounds;
use crate::feature::record_buf::Strand;

#[derive(Clone, Eq, PartialEq)]
pub(crate) struct Fields<const N: usize> {
    pub(crate) buf: Vec<u8>,
    pub(crate) bounds: Bounds<N>,
}

impl<const N: usize> Fields<N> {
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

    pub(super) fn feature_end(&self) -> Option<io::Result<Position>> {
        const MISSING: &[u8] = b"0";

        let src = &self.buf[self.bounds.feature_end_range()];

        match src {
            MISSING => None,
            _ => Some(parse_int::<usize>(src).and_then(|n| {
                Position::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })),
        }
    }

    pub(super) fn name(&self) -> Option<&[u8]> {
        const NAME_INDEX: usize = 0;
        self.get(NAME_INDEX)
    }

    pub(super) fn score(&self) -> Option<io::Result<u16>> {
        const SCORE_INDEX: usize = 1;
        self.get(SCORE_INDEX).map(parse_int)
    }

    pub(super) fn strand(&self) -> Option<io::Result<Option<Strand>>> {
        const STRAND_INDEX: usize = 2;
        self.get(STRAND_INDEX).map(parse_strand)
    }

    pub(super) fn get(&self, i: usize) -> Option<&[u8]> {
        self.bounds.get(i).map(|range| &self.buf[range])
    }
}

impl Default for Fields<3> {
    fn default() -> Self {
        Self {
            buf: Vec::from(*b"sq001"),
            bounds: Bounds::default(),
        }
    }
}

impl Default for Fields<4> {
    fn default() -> Self {
        Self {
            buf: Vec::from(*b"sq001."),
            bounds: Bounds::default(),
        }
    }
}

impl Default for Fields<5> {
    fn default() -> Self {
        Self {
            buf: Vec::from(*b"sq001.0"),
            bounds: Bounds::default(),
        }
    }
}

fn parse_int<N: FromLexical>(buf: &[u8]) -> io::Result<N> {
    lexical_core::parse(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_strand(buf: &[u8]) -> io::Result<Option<Strand>> {
    const MISSING: &[u8] = b".";
    const FORWARD: &[u8] = b"+";
    const REVERSE: &[u8] = b"-";

    match buf {
        MISSING => Ok(None),
        FORWARD => Ok(Some(Strand::Forward)),
        REVERSE => Ok(Some(Strand::Reverse)),
        _ => Err(io::Error::new(io::ErrorKind::InvalidData, "invalid strand")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_strand() -> io::Result<()> {
        assert!(parse_strand(b".")?.is_none());
        assert_eq!(parse_strand(b"+")?, Some(Strand::Forward));
        assert_eq!(parse_strand(b"-")?, Some(Strand::Reverse));

        assert!(matches!(parse_strand(b"n"), Err(e) if e.kind() == io::ErrorKind::InvalidData));

        Ok(())
    }
}
