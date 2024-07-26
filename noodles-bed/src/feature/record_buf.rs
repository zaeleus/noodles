//! BED record and fields.

pub mod builder;
pub mod color;
pub mod strand;

pub use self::{builder::Builder, color::Color, strand::Strand};

use bstr::{BStr, BString};
use noodles_core::Position;

#[derive(Clone, Debug, Eq, PartialEq)]
struct StandardFields<const N: usize> {
    reference_sequence_name: BString,
    feature_start: Position,
    feature_end: Option<Position>,
    name: Option<BString>,
    score: u16,
    strand: Option<Strand>,
}

impl<const N: usize> StandardFields<N> {
    const fn len(&self) -> usize {
        N
    }
}

impl<const N: usize> Default for StandardFields<N> {
    fn default() -> Self {
        Self {
            reference_sequence_name: BString::default(),
            feature_start: Position::MIN,
            feature_end: None,
            name: None,
            score: 0,
            strand: None,
        }
    }
}

/// A feature record buffer.
#[derive(Clone, Default, Debug, Eq, PartialEq)]
pub struct RecordBuf<const N: usize> {
    standard_fields: StandardFields<N>,
    other_fields: Vec<Option<Vec<u8>>>,
}

impl<const N: usize> RecordBuf<N> {
    /// Returns a feature record builder.
    pub fn builder() -> Builder<N> {
        Builder::default()
    }

    /// Returns the number of standard fields.
    pub const fn standard_field_count(&self) -> usize {
        self.standard_fields.len()
    }

    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &[u8] {
        &self.standard_fields.reference_sequence_name
    }

    /// Returns the feature start.
    pub fn feature_start(&self) -> Position {
        self.standard_fields.feature_start
    }

    /// Returns the feature end.
    pub fn feature_end(&self) -> Option<Position> {
        self.standard_fields.feature_end
    }
}

impl RecordBuf<4> {
    /// Returns the name.
    pub fn name(&self) -> Option<&BStr> {
        self.standard_fields.name.as_ref().map(|name| name.as_ref())
    }
}

impl RecordBuf<5> {
    /// Returns the name.
    pub fn name(&self) -> Option<&BStr> {
        self.standard_fields.name.as_ref().map(|name| name.as_ref())
    }

    /// Returns the score.
    pub fn score(&self) -> u16 {
        self.standard_fields.score
    }
}

impl RecordBuf<6> {
    /// Returns the name.
    pub fn name(&self) -> Option<&BStr> {
        self.standard_fields.name.as_ref().map(|name| name.as_ref())
    }

    /// Returns the score.
    pub fn score(&self) -> u16 {
        self.standard_fields.score
    }

    /// Returns the strand.
    pub fn strand(&self) -> Option<Strand> {
        self.standard_fields.strand
    }
}
