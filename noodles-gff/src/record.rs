//! GFF record.

pub mod attributes;
pub(crate) mod fields;

use std::{fmt, io};

use bstr::{BStr, ByteSlice};
use noodles_core::Position;

pub use self::attributes::Attributes;
use self::fields::Fields;
use crate::feature::record::{Phase, Strand};

const MISSING: &[u8] = b".";

/// A GFF record.
#[derive(Clone, Eq, PartialEq)]
pub struct Record<'l>(Fields<'l>);

impl<'l> Record<'l> {
    pub(super) fn try_new(src: &'l [u8]) -> io::Result<Self> {
        Fields::try_new(src).map(Self)
    }

    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &BStr {
        self.0.reference_sequence_name()
    }

    /// Returns the source.
    pub fn source(&self) -> &BStr {
        self.0.source()
    }

    /// Returns the feature type.
    pub fn ty(&self) -> &BStr {
        self.0.ty()
    }

    /// Returns the start position.
    pub fn start(&self) -> io::Result<Position> {
        self.0.start()
    }

    /// Returns the end position.
    pub fn end(&self) -> io::Result<Position> {
        self.0.end()
    }

    /// Returns the score.
    pub fn score(&self) -> Option<io::Result<f32>> {
        parse_score(self.0.score())
    }

    /// Returns the strand.
    pub fn strand(&self) -> io::Result<Strand> {
        parse_strand(self.0.strand())
    }

    /// Returns the phase.
    pub fn phase(&self) -> Option<io::Result<Phase>> {
        parse_phase(self.0.phase())
    }

    /// Returns the attributes.
    pub fn attributes(&self) -> Attributes<'_> {
        self.0.attributes()
    }
}

impl fmt::Debug for Record<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Record")
            .field("reference_sequence_name", &self.reference_sequence_name())
            .field("source", &self.source())
            .field("ty", &self.ty())
            .field("start", &self.start())
            .field("end", &self.end())
            .field("score", &self.score())
            .field("strand", &self.strand())
            .field("phase", &self.phase())
            .field("attributes", &self.attributes())
            .finish()
    }
}

impl super::feature::Record for Record<'_> {
    fn reference_sequence_name(&self) -> &BStr {
        self.reference_sequence_name().as_bytes().as_bstr()
    }

    fn source(&self) -> &BStr {
        self.source().as_bytes().as_bstr()
    }

    fn ty(&self) -> &BStr {
        self.ty().as_bytes().as_bstr()
    }

    fn feature_start(&self) -> io::Result<Position> {
        self.start()
    }

    fn feature_end(&self) -> io::Result<Position> {
        self.end()
    }

    fn score(&self) -> Option<io::Result<f32>> {
        self.score()
    }

    fn strand(&self) -> io::Result<Strand> {
        self.strand()
    }

    fn phase(&self) -> Option<io::Result<Phase>> {
        self.phase()
    }

    fn attributes(&self) -> Box<dyn super::feature::record::Attributes + '_> {
        Box::new(self.attributes())
    }
}

fn parse_score(src: &[u8]) -> Option<io::Result<f32>> {
    match src {
        MISSING => None,
        _ => Some(
            lexical_core::parse(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        ),
    }
}

fn parse_strand(src: &[u8]) -> io::Result<Strand> {
    match src {
        b"." => Ok(Strand::None),
        b"+" => Ok(Strand::Forward),
        b"-" => Ok(Strand::Reverse),
        b"?" => Ok(Strand::Unknown),
        _ => Err(io::Error::new(io::ErrorKind::InvalidData, "invalid strand")),
    }
}

fn parse_phase(src: &[u8]) -> Option<io::Result<Phase>> {
    match src {
        MISSING => None,
        b"0" => Some(Ok(Phase::Zero)),
        b"1" => Some(Ok(Phase::One)),
        b"2" => Some(Ok(Phase::Two)),
        _ => Some(Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid phase",
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_score() -> io::Result<()> {
        assert!(parse_score(b".").is_none());
        assert_eq!(parse_score(b"0.0").transpose()?, Some(0.0));

        assert!(matches!(
            parse_phase(b""),
            Some(Err(e)) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_parse_strand() -> io::Result<()> {
        assert_eq!(parse_strand(b".")?, Strand::None);
        assert_eq!(parse_strand(b"+")?, Strand::Forward);
        assert_eq!(parse_strand(b"-")?, Strand::Reverse);
        assert_eq!(parse_strand(b"?")?, Strand::Unknown);

        assert!(matches!(
            parse_strand(b""),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_parse_phase() -> io::Result<()> {
        assert!(parse_phase(b".").is_none());
        assert_eq!(parse_phase(b"0").transpose()?, Some(Phase::Zero));
        assert_eq!(parse_phase(b"1").transpose()?, Some(Phase::One));
        assert_eq!(parse_phase(b"2").transpose()?, Some(Phase::Two));

        assert!(matches!(
            parse_phase(b""),
            Some(Err(e)) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
