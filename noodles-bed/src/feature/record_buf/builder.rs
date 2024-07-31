//! BED record builder.

use bstr::BString;
use noodles_core::Position;

use super::{OtherFields, RecordBuf, StandardFields, Strand};

/// A feature record builder.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Builder<const N: usize> {
    reference_sequence_name: BString,
    feature_start: Option<Position>,
    feature_end: Option<Position>,
    name: Option<BString>,
    score: u16,
    strand: Option<Strand>,
    other_fields: OtherFields,
}

impl<const N: usize> Builder<N> {
    /// Builds a feature record buffer.
    pub fn build(self) -> RecordBuf<N> {
        RecordBuf {
            standard_fields: StandardFields {
                reference_sequence_name: self.reference_sequence_name,
                feature_start: self.feature_start.unwrap_or(Position::MIN),
                feature_end: self.feature_end,
                name: self.name,
                score: self.score,
                strand: self.strand,
            },
            other_fields: self.other_fields,
        }
    }
}

impl Builder<3> {
    /// Sets the reference sequence name.
    pub fn set_reference_sequence_name<M>(mut self, reference_sequence_name: M) -> Self
    where
        M: Into<BString>,
    {
        self.reference_sequence_name = reference_sequence_name.into();
        self
    }

    /// Sets the feature start position.
    pub fn set_feature_start(mut self, start_position: Position) -> Self {
        self.feature_start = Some(start_position);
        self
    }

    /// Sets the feature end position.
    pub fn set_feature_end(mut self, end_position: Position) -> Self {
        self.feature_end = Some(end_position);
        self
    }

    /// Sets the list of raw optional fields.
    pub fn set_other_fields(mut self, other_fields: OtherFields) -> Self {
        self.other_fields = other_fields;
        self
    }
}

impl Builder<4> {
    /// Sets the reference sequence name.
    pub fn set_reference_sequence_name<M>(mut self, reference_sequence_name: M) -> Self
    where
        M: Into<BString>,
    {
        self.reference_sequence_name = reference_sequence_name.into();
        self
    }

    /// Sets the feature start position.
    pub fn set_feature_start(mut self, start_position: Position) -> Self {
        self.feature_start = Some(start_position);
        self
    }

    /// Sets the feature end position.
    pub fn set_feature_end(mut self, end_position: Position) -> Self {
        self.feature_end = Some(end_position);
        self
    }

    /// Sets the name.
    pub fn set_name<M>(mut self, name: M) -> Self
    where
        M: Into<BString>,
    {
        self.name = Some(name.into());
        self
    }

    /// Sets the list of raw optional fields.
    pub fn set_other_fields(mut self, other_fields: OtherFields) -> Self {
        self.other_fields = other_fields;
        self
    }
}

impl Builder<5> {
    /// Sets the reference sequence name.
    pub fn set_reference_sequence_name<M>(mut self, reference_sequence_name: M) -> Self
    where
        M: Into<BString>,
    {
        self.reference_sequence_name = reference_sequence_name.into();
        self
    }

    /// Sets the feature start position.
    pub fn set_feature_start(mut self, start_position: Position) -> Self {
        self.feature_start = Some(start_position);
        self
    }

    /// Sets the feature end position.
    pub fn set_feature_end(mut self, end_position: Position) -> Self {
        self.feature_end = Some(end_position);
        self
    }

    /// Sets the name.
    pub fn set_name<M>(mut self, name: M) -> Self
    where
        M: Into<BString>,
    {
        self.name = Some(name.into());
        self
    }

    /// Sets the score.
    pub fn set_score(mut self, score: u16) -> Self {
        self.score = score;
        self
    }

    /// Sets the list of raw optional fields.
    pub fn set_other_fields(mut self, other_fields: OtherFields) -> Self {
        self.other_fields = other_fields;
        self
    }
}

impl Builder<6> {
    /// Sets the reference sequence name.
    pub fn set_reference_sequence_name<M>(mut self, reference_sequence_name: M) -> Self
    where
        M: Into<BString>,
    {
        self.reference_sequence_name = reference_sequence_name.into();
        self
    }

    /// Sets the feature start position.
    pub fn set_feature_start(mut self, start_position: Position) -> Self {
        self.feature_start = Some(start_position);
        self
    }

    /// Sets the feature end position.
    pub fn set_feature_end(mut self, end_position: Position) -> Self {
        self.feature_end = Some(end_position);
        self
    }

    /// Sets the name.
    pub fn set_name<M>(mut self, name: M) -> Self
    where
        M: Into<BString>,
    {
        self.name = Some(name.into());
        self
    }

    /// Sets the score.
    pub fn set_score(mut self, score: u16) -> Self {
        self.score = score;
        self
    }

    /// Sets the strand.
    pub fn set_strand(mut self, strand: Strand) -> Self {
        self.strand = Some(strand);
        self
    }

    /// Sets the list of raw optional fields.
    pub fn set_other_fields(mut self, other_fields: OtherFields) -> Self {
        self.other_fields = other_fields;
        self
    }
}
