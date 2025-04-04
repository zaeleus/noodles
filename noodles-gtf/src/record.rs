//! GTF record.

mod attributes;
mod fields;

use std::{fmt, io};

use bstr::{BStr, ByteSlice};
use noodles_core::Position;
use noodles_gff::{
    self as gff,
    feature::record::{Phase, Strand},
};

pub use self::attributes::Attributes;
use self::fields::Fields;

/// A GTF record.
#[derive(Clone, Eq, PartialEq)]
pub struct Record<'l>(Fields<'l>);

const MISSING: &str = ".";

impl<'l> Record<'l> {
    pub(super) fn try_new(src: &'l str) -> io::Result<Self> {
        Fields::try_new(src).map(Self)
    }

    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &str {
        self.0.reference_sequence_name()
    }

    /// Returns the source.
    pub fn source(&self) -> &str {
        self.0.source()
    }

    /// Returns the feature type.
    pub fn ty(&self) -> &str {
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
    #[deprecated(since = "0.41.0", note = "Use `Record::phase` instead.")]
    pub fn frame(&self) -> Option<io::Result<Phase>> {
        self.phase()
    }

    /// Returns the phase.
    pub fn phase(&self) -> Option<io::Result<Phase>> {
        parse_phase(self.0.phase())
    }

    /// Returns the attributes.
    pub fn attributes(&self) -> io::Result<Attributes<'_>> {
        Attributes::try_new(self.0.attributes().as_bytes())
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

impl gff::feature::Record for Record<'_> {
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

    fn attributes(&self) -> Box<dyn gff::feature::record::Attributes + '_> {
        Box::new(self.attributes().unwrap()) // TODO
    }
}

fn parse_score(s: &str) -> Option<io::Result<f32>> {
    match s {
        MISSING => None,
        _ => Some(
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        ),
    }
}

fn parse_strand(s: &str) -> io::Result<Strand> {
    match s {
        MISSING => Ok(Strand::None),
        "+" => Ok(Strand::Forward),
        "-" => Ok(Strand::Reverse),
        _ => Err(io::Error::new(io::ErrorKind::InvalidData, "invalid strand")),
    }
}

fn parse_phase(s: &str) -> Option<io::Result<Phase>> {
    match s {
        MISSING => None,
        "0" => Some(Ok(Phase::Zero)),
        "1" => Some(Ok(Phase::One)),
        "2" => Some(Ok(Phase::Two)),
        _ => Some(Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid phase",
        ))),
    }
}
