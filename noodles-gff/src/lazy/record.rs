mod attributes;
mod bounds;

use self::attributes::Attributes;
pub(crate) use self::bounds::Bounds;

const MISSING: &str = ".";

/// An immutable, lazily-evalulated GFF record.
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
    pub fn start(&self) -> &str {
        &self.buf[self.bounds.start_range()]
    }

    /// Returns the end position.
    pub fn end(&self) -> &str {
        &self.buf[self.bounds.end_range()]
    }

    /// Returns the score.
    pub fn score(&self) -> &str {
        &self.buf[self.bounds.score_range()]
    }

    /// Returns the strand.
    pub fn strand(&self) -> &str {
        &self.buf[self.bounds.strand_range()]
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

impl From<Record> for String {
    fn from(record: Record) -> Self {
        record.buf
    }
}
