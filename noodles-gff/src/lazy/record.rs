mod attributes;
mod bounds;
mod position;
mod strand;

use core::fmt;

pub(crate) use self::bounds::Bounds;
use self::{attributes::Attributes, position::Position, strand::Strand};

const MISSING: &str = ".";

/// An immutable, lazily-evalulated GFF record.
#[derive(Clone, Eq, PartialEq)]
pub struct Record {
    pub(crate) buf: String,
    pub(crate) bounds: Bounds,
}

impl Record {
    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &str {
        &self.buf[self.bounds.reference_sequence_name_range()]
    }

    /// Returns the source.
    pub fn source(&self) -> &str {
        &self.buf[self.bounds.source_range()]
    }

    /// Returns the feature type.
    pub fn ty(&self) -> &str {
        &self.buf[self.bounds.type_range()]
    }

    /// Returns the start position.
    pub fn start(&self) -> Position<'_> {
        let buf = &self.buf[self.bounds.start_range()];
        Position::new(buf)
    }

    /// Returns the end position.
    pub fn end(&self) -> Position<'_> {
        let buf = &self.buf[self.bounds.end_range()];
        Position::new(buf)
    }

    /// Returns the score.
    pub fn score(&self) -> &str {
        &self.buf[self.bounds.score_range()]
    }

    /// Returns the strand.
    pub fn strand(&self) -> Strand<'_> {
        let buf = &self.buf[self.bounds.strand_range()];
        Strand::new(buf)
    }

    /// Returns the phase.
    pub fn phase(&self) -> &str {
        &self.buf[self.bounds.phase_range()]
    }

    /// Returns the attributes.
    pub fn attributes(&self) -> Attributes<'_> {
        match &self.buf[self.bounds.attributes_range()] {
            MISSING => Attributes::new(""),
            buf => Attributes::new(buf),
        }
    }
}

impl fmt::Debug for Record {
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

impl Default for Record {
    fn default() -> Self {
        Self {
            buf: String::from("...11...."),
            bounds: Bounds::default(),
        }
    }
}

impl From<Record> for String {
    fn from(record: Record) -> Self {
        record.buf
    }
}
