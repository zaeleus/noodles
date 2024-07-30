//! BED record and fields.

pub mod builder;
pub mod color;
mod other_fields;
pub mod strand;

use std::io;

use bstr::{BStr, BString};
use noodles_core::Position;

pub use self::{builder::Builder, color::Color, other_fields::OtherFields, strand::Strand};

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
    other_fields: OtherFields,
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

    /// Returns the other fields.
    pub fn other_fields(&self) -> &OtherFields {
        &self.other_fields
    }
}

impl RecordBuf<3> {
    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &BStr {
        self.standard_fields.reference_sequence_name.as_ref()
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

impl super::Record<3> for RecordBuf<3> {
    fn reference_sequence_name(&self) -> &BStr {
        self.reference_sequence_name()
    }

    fn feature_start(&self) -> io::Result<Position> {
        Ok(self.feature_start())
    }

    fn feature_end(&self) -> Option<io::Result<Position>> {
        self.feature_end().map(Ok)
    }

    fn name(&self) -> Option<Option<&BStr>> {
        None
    }

    fn score(&self) -> Option<io::Result<u16>> {
        None
    }

    fn strand(&self) -> Option<io::Result<Option<Strand>>> {
        None
    }

    fn other_fields(&self) -> Box<dyn super::record::OtherFields + '_> {
        Box::new(self.other_fields())
    }
}

impl RecordBuf<4> {
    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &BStr {
        self.standard_fields.reference_sequence_name.as_ref()
    }

    /// Returns the feature start.
    pub fn feature_start(&self) -> Position {
        self.standard_fields.feature_start
    }

    /// Returns the feature end.
    pub fn feature_end(&self) -> Option<Position> {
        self.standard_fields.feature_end
    }

    /// Returns the name.
    pub fn name(&self) -> Option<&BStr> {
        self.standard_fields.name.as_ref().map(|name| name.as_ref())
    }
}

impl super::Record<4> for RecordBuf<4> {
    fn reference_sequence_name(&self) -> &BStr {
        self.reference_sequence_name()
    }

    fn feature_start(&self) -> io::Result<Position> {
        Ok(self.feature_start())
    }

    fn feature_end(&self) -> Option<io::Result<Position>> {
        self.feature_end().map(Ok)
    }

    fn name(&self) -> Option<Option<&BStr>> {
        Some(self.name())
    }

    fn score(&self) -> Option<io::Result<u16>> {
        None
    }

    fn strand(&self) -> Option<io::Result<Option<Strand>>> {
        None
    }

    fn other_fields(&self) -> Box<dyn super::record::OtherFields + '_> {
        Box::new(self.other_fields())
    }
}

impl RecordBuf<5> {
    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &BStr {
        self.standard_fields.reference_sequence_name.as_ref()
    }

    /// Returns the feature start.
    pub fn feature_start(&self) -> Position {
        self.standard_fields.feature_start
    }

    /// Returns the feature end.
    pub fn feature_end(&self) -> Option<Position> {
        self.standard_fields.feature_end
    }

    /// Returns the name.
    pub fn name(&self) -> Option<&BStr> {
        self.standard_fields.name.as_ref().map(|name| name.as_ref())
    }

    /// Returns the score.
    pub fn score(&self) -> u16 {
        self.standard_fields.score
    }
}

impl super::Record<5> for RecordBuf<5> {
    fn reference_sequence_name(&self) -> &BStr {
        self.reference_sequence_name()
    }

    fn feature_start(&self) -> io::Result<Position> {
        Ok(self.feature_start())
    }

    fn feature_end(&self) -> Option<io::Result<Position>> {
        self.feature_end().map(Ok)
    }

    fn name(&self) -> Option<Option<&BStr>> {
        Some(self.name())
    }

    fn score(&self) -> Option<io::Result<u16>> {
        Some(Ok(self.score()))
    }

    fn strand(&self) -> Option<io::Result<Option<Strand>>> {
        None
    }

    fn other_fields(&self) -> Box<dyn super::record::OtherFields + '_> {
        Box::new(self.other_fields())
    }
}

impl RecordBuf<6> {
    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &BStr {
        self.standard_fields.reference_sequence_name.as_ref()
    }

    /// Returns the feature start.
    pub fn feature_start(&self) -> Position {
        self.standard_fields.feature_start
    }

    /// Returns the feature end.
    pub fn feature_end(&self) -> Option<Position> {
        self.standard_fields.feature_end
    }

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

impl super::Record<6> for RecordBuf<6> {
    fn reference_sequence_name(&self) -> &BStr {
        self.reference_sequence_name()
    }

    fn feature_start(&self) -> io::Result<Position> {
        Ok(self.feature_start())
    }

    fn feature_end(&self) -> Option<io::Result<Position>> {
        self.feature_end().map(Ok)
    }

    fn name(&self) -> Option<Option<&BStr>> {
        Some(self.name())
    }

    fn score(&self) -> Option<io::Result<u16>> {
        Some(Ok(self.score()))
    }

    fn strand(&self) -> Option<io::Result<Option<Strand>>> {
        Some(Ok(self.strand()))
    }

    fn other_fields(&self) -> Box<dyn super::record::OtherFields + '_> {
        Box::new(self.other_fields())
    }
}
