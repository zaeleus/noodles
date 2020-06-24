mod attributes;
mod strand;

pub use self::{attributes::Attributes, strand::Strand};

use std::{error, fmt};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Header {
    SeqName,
    Source,
    Feature,
    Start,
    End,
    Score,
    Strand,
    Frame,
    Attributes,
}

#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    Missing(Header),
    Empty(Header),
    Parse(Header, String),
}

impl error::Error for Error {}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

pub struct Record {
    reference_sequence_name: String,
    source: String,
    feature: String,
    start: String,
    end: String,
    score: f32,
    strand: String,
    frame: String,
    attributes: Option<String>,
}

impl Record {
    pub fn reference_sequence_name(&self) -> &str {
        &self.reference_sequence_name
    }

    pub fn source(&self) -> &str {
        &self.source
    }

    pub fn feature(&self) -> &str {
        &self.feature
    }

    pub fn start(&self) -> &str {
        &self.start
    }

    pub fn end(&self) -> &str {
        &self.end
    }

    pub fn score(&self) -> f32 {
        self.score
    }

    pub fn strand(&self) -> &str {
        &self.strand
    }

    pub fn frame(&self) -> &str {
        &self.frame
    }

    pub fn attributes(&self) -> Option<&str> {
        self.attributes.as_deref()
    }
}
