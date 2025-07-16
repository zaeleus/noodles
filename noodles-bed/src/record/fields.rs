mod bounds;

use std::io;

use bstr::{BStr, ByteSlice};
use lexical_core::FromLexical;
use noodles_core::Position;

pub(crate) use self::bounds::Bounds;
use crate::feature::record::Strand;

#[derive(Clone, Eq, PartialEq)]
pub(crate) struct Fields<const N: usize> {
    pub(crate) buf: Vec<u8>,
    pub(crate) bounds: Bounds<N>,
}

impl<const N: usize> Fields<N> {
    pub(super) fn get(&self, i: usize) -> Option<&[u8]> {
        self.bounds.get(i).map(|range| &self.buf[range])
    }
}

impl Fields<3> {
    pub(super) fn reference_sequence_name(&self) -> &BStr {
        let src = &self.buf[self.bounds.reference_sequence_name_range()];
        parse_reference_sequence_name(src)
    }

    pub(super) fn feature_start(&self) -> io::Result<Position> {
        let src = &self.buf[self.bounds.feature_start_range()];
        parse_feature_start(src)
    }

    pub(super) fn feature_end(&self) -> Option<io::Result<Position>> {
        let src = &self.buf[self.bounds.feature_end_range()];
        parse_feature_end(src)
    }
}

impl Fields<4> {
    pub(super) fn reference_sequence_name(&self) -> &BStr {
        let src = &self.buf[self.bounds.reference_sequence_name_range()];
        parse_reference_sequence_name(src)
    }

    pub(super) fn feature_start(&self) -> io::Result<Position> {
        let src = &self.buf[self.bounds.feature_start_range()];
        parse_feature_start(src)
    }

    pub(super) fn feature_end(&self) -> Option<io::Result<Position>> {
        let src = &self.buf[self.bounds.feature_end_range()];
        parse_feature_end(src)
    }

    pub(super) fn name(&self) -> Option<&BStr> {
        let src = &self.buf[self.bounds.name_range()];
        parse_name(src)
    }
}

impl Fields<5> {
    pub(super) fn reference_sequence_name(&self) -> &BStr {
        let src = &self.buf[self.bounds.reference_sequence_name_range()];
        parse_reference_sequence_name(src)
    }

    pub(super) fn feature_start(&self) -> io::Result<Position> {
        let src = &self.buf[self.bounds.feature_start_range()];
        parse_feature_start(src)
    }

    pub(super) fn feature_end(&self) -> Option<io::Result<Position>> {
        let src = &self.buf[self.bounds.feature_end_range()];
        parse_feature_end(src)
    }

    pub(super) fn name(&self) -> Option<&BStr> {
        let src = &self.buf[self.bounds.name_range()];
        parse_name(src)
    }

    pub(super) fn score(&self) -> io::Result<u16> {
        let src = &self.buf[self.bounds.score_range()];
        parse_int(src)
    }
}

impl Fields<6> {
    pub(super) fn reference_sequence_name(&self) -> &BStr {
        let src = &self.buf[self.bounds.reference_sequence_name_range()];
        parse_reference_sequence_name(src)
    }

    pub(super) fn feature_start(&self) -> io::Result<Position> {
        let src = &self.buf[self.bounds.feature_start_range()];
        parse_feature_start(src)
    }

    pub(super) fn feature_end(&self) -> Option<io::Result<Position>> {
        let src = &self.buf[self.bounds.feature_end_range()];
        parse_feature_end(src)
    }

    pub(super) fn name(&self) -> Option<&BStr> {
        let src = &self.buf[self.bounds.name_range()];
        parse_name(src)
    }

    pub(super) fn score(&self) -> io::Result<u16> {
        let src = &self.buf[self.bounds.score_range()];
        parse_int(src)
    }

    pub(super) fn strand(&self) -> io::Result<Option<Strand>> {
        let src = &self.buf[self.bounds.strand_range()];
        parse_strand(src)
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

impl Default for Fields<6> {
    fn default() -> Self {
        Self {
            buf: Vec::from(*b"sq001.0."),
            bounds: Bounds::default(),
        }
    }
}

fn parse_int<N: FromLexical>(buf: &[u8]) -> io::Result<N> {
    lexical_core::parse(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_reference_sequence_name(buf: &[u8]) -> &BStr {
    buf.as_bstr()
}

fn parse_feature_start(buf: &[u8]) -> io::Result<Position> {
    parse_int::<usize>(buf).and_then(|n| {
        n.checked_add(1)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "attempt to add with overflow")
            })
            .and_then(|m| {
                Position::try_from(m).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    })
}

fn parse_feature_end(buf: &[u8]) -> Option<io::Result<Position>> {
    const MISSING: &[u8] = b"0";

    match buf {
        MISSING => None,
        _ => Some(parse_int::<usize>(buf).and_then(|n| {
            Position::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })),
    }
}

fn parse_name(buf: &[u8]) -> Option<&BStr> {
    const MISSING: &[u8] = b".";

    match buf {
        MISSING => None,
        _ => Some(buf.as_bstr()),
    }
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
    fn test_parse_feature_end() -> io::Result<()> {
        assert!(parse_feature_end(b"0").is_none());
        assert_eq!(parse_feature_end(b"1").transpose()?, Some(Position::MIN));

        assert_eq!(
            parse_feature_end(Position::MAX.to_string().as_bytes()).transpose()?,
            Some(Position::MAX)
        );

        assert!(matches!(
            parse_feature_end(b"-1").transpose(),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_parse_name() {
        assert!(parse_name(b".").is_none());
        assert_eq!(parse_name(b"ndls"), Some(BStr::new("ndls")));
    }

    #[test]
    fn test_parse_strand() -> io::Result<()> {
        assert!(parse_strand(b".")?.is_none());
        assert_eq!(parse_strand(b"+")?, Some(Strand::Forward));
        assert_eq!(parse_strand(b"-")?, Some(Strand::Reverse));

        assert!(matches!(parse_strand(b"n"), Err(e) if e.kind() == io::ErrorKind::InvalidData));

        Ok(())
    }
}
