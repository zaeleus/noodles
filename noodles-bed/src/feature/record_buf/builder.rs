//! BED record builder.

use bstr::BString;
use noodles_core::Position;

use super::{OtherFields, RecordBuf, StandardFields, Strand};

/// A feature record builder.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Builder<const N: usize> {
    reference_sequence_name: BString,
    feature_start: Option<Position>,
    feature_end: Option<Position>,
    name: Option<Vec<u8>>,
    score: Option<u16>,
    strand: Option<Strand>,
    other_fields: OtherFields,
}

impl Builder<3> {
    /// Sets the reference sequence name (`chrom`).
    pub fn set_reference_sequence_name<M>(mut self, reference_sequence_name: M) -> Self
    where
        M: Into<BString>,
    {
        self.reference_sequence_name = reference_sequence_name.into();
        self
    }

    /// Sets the feature start position (`chromStart`).
    pub fn set_start_position(mut self, start_position: Position) -> Self {
        self.feature_start = Some(start_position);
        self
    }

    /// Sets the feature end position (`chromEnd`).
    pub fn set_end_position(mut self, end_position: Position) -> Self {
        self.feature_end = Some(end_position);
        self
    }

    /// Sets the list of raw optional fields.
    pub fn set_other_fields(mut self, other_fields: OtherFields) -> Self {
        self.other_fields = other_fields;
        self
    }
}

impl Builder<3> {
    /// Builds a BED3 record.
    pub fn build(self) -> RecordBuf<3> {
        RecordBuf {
            standard_fields: StandardFields {
                reference_sequence_name: self.reference_sequence_name,
                feature_start: self.feature_start.unwrap_or(Position::MIN),
                feature_end: self.feature_end,
                name: None,
                score: 0,
                strand: None,
            },
            other_fields: self.other_fields,
        }
    }
}
