mod bounds;

use std::io;

use noodles_core::Position;

pub(crate) use self::bounds::Bounds;
use super::Attributes;
use crate::record::Strand;

#[derive(Clone, Eq, PartialEq)]
pub(crate) struct Fields {
    pub(crate) buf: String,
    pub(crate) bounds: Bounds,
}

impl Fields {
    pub fn reference_sequence_name(&self) -> &str {
        &self.buf[self.bounds.reference_sequence_name_range()]
    }

    pub fn source(&self) -> &str {
        &self.buf[self.bounds.source_range()]
    }

    pub fn ty(&self) -> &str {
        &self.buf[self.bounds.type_range()]
    }

    pub fn start(&self) -> io::Result<Position> {
        let src = &self.buf[self.bounds.start_range()];
        parse_position(src)
    }

    pub fn end(&self) -> io::Result<Position> {
        let src = &self.buf[self.bounds.end_range()];
        parse_position(src)
    }

    pub fn score(&self) -> &str {
        &self.buf[self.bounds.score_range()]
    }

    pub fn strand(&self) -> io::Result<Strand> {
        let src = &self.buf[self.bounds.strand_range()];
        parse_strand(src)
    }

    pub fn phase(&self) -> &str {
        &self.buf[self.bounds.phase_range()]
    }

    pub fn attributes(&self) -> Attributes<'_> {
        const MISSING: &str = ".";

        match &self.buf[self.bounds.attributes_range()] {
            MISSING => Attributes::new(""),
            buf => Attributes::new(buf),
        }
    }
}

impl Default for Fields {
    fn default() -> Self {
        Self {
            buf: String::from("...11...."),
            bounds: Bounds::default(),
        }
    }
}

fn parse_position(s: &str) -> io::Result<Position> {
    s.parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_strand(s: &str) -> io::Result<Strand> {
    s.parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
