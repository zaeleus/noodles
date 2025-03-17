//! GTF record.

mod attributes;
mod fields;

use std::io;

use noodles_core::Position;
use noodles_gff::feature::record::Strand;

pub use self::attributes::Attributes;
use self::fields::Fields;
use super::record_buf::Frame;

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

    /// Returns the frame.
    pub fn frame(&self) -> Option<io::Result<Frame>> {
        parse_frame(self.0.frame())
    }

    /// Returns the attributes.
    pub fn attributes(&self) -> io::Result<Attributes<'_>> {
        Attributes::try_new(self.0.attributes())
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

fn parse_frame(s: &str) -> Option<io::Result<Frame>> {
    match s {
        MISSING => None,
        _ => Some(
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        ),
    }
}
