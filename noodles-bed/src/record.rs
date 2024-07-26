pub(crate) mod fields;
mod other_fields;

use std::{fmt, io};

use bstr::BStr;
use noodles_core::Position;

use self::fields::Fields;
pub use self::other_fields::OtherFields;
use crate::feature::record_buf::Strand;

/// A BED record.
#[derive(Clone, Eq, PartialEq)]
pub struct Record<const N: usize>(pub(crate) Fields<N>);

impl<const N: usize> Record<N> {
    /// Returns the number of standard fields.
    pub const fn standard_field_count(&self) -> usize {
        N
    }

    /// Returns the other fields.
    pub fn other_fields(&self) -> OtherFields<'_, N> {
        OtherFields::new(&self.0)
    }

    /// Returns the number of fields.
    ///
    /// This is guaranteed to be >= 3.
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        const MIN_FIELD_COUNT: usize = 3;
        MIN_FIELD_COUNT + self.other_fields().len()
    }
}

impl Record<3> {
    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &BStr {
        self.0.reference_sequence_name()
    }

    /// Returns the feature start.
    pub fn feature_start(&self) -> io::Result<Position> {
        self.0.feature_start()
    }

    /// Returns the feature end.
    pub fn feature_end(&self) -> Option<io::Result<Position>> {
        self.0.feature_end()
    }
}

impl Record<4> {
    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &BStr {
        self.0.reference_sequence_name()
    }

    /// Returns the feature start.
    pub fn feature_start(&self) -> io::Result<Position> {
        self.0.feature_start()
    }

    /// Returns the feature end.
    pub fn feature_end(&self) -> Option<io::Result<Position>> {
        self.0.feature_end()
    }

    /// Returns the name.
    pub fn name(&self) -> Option<&BStr> {
        self.0.name()
    }
}

impl Record<5> {
    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &BStr {
        self.0.reference_sequence_name()
    }

    /// Returns the feature start.
    pub fn feature_start(&self) -> io::Result<Position> {
        self.0.feature_start()
    }

    /// Returns the feature end.
    pub fn feature_end(&self) -> Option<io::Result<Position>> {
        self.0.feature_end()
    }

    /// Returns the name.
    pub fn name(&self) -> Option<&BStr> {
        self.0.name()
    }

    /// Returns the score.
    pub fn score(&self) -> io::Result<u16> {
        self.0.score()
    }
}

impl Record<6> {
    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &BStr {
        self.0.reference_sequence_name()
    }

    /// Returns the feature start.
    pub fn feature_start(&self) -> io::Result<Position> {
        self.0.feature_start()
    }

    /// Returns the feature end.
    pub fn feature_end(&self) -> Option<io::Result<Position>> {
        self.0.feature_end()
    }

    /// Returns the name.
    pub fn name(&self) -> Option<&BStr> {
        self.0.name()
    }

    /// Returns the score.
    pub fn score(&self) -> io::Result<u16> {
        self.0.score()
    }

    /// Returns the strand.
    pub fn strand(&self) -> io::Result<Option<Strand>> {
        self.0.strand()
    }
}

impl fmt::Debug for Record<3> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Record")
            .field("reference_sequence_name", &self.reference_sequence_name())
            .field("feature_start", &self.feature_start())
            .field("feature_end", &self.feature_end())
            .finish_non_exhaustive()
    }
}

impl Default for Record<3> {
    fn default() -> Self {
        Self(Fields::default())
    }
}

impl Default for Record<4> {
    fn default() -> Self {
        Self(Fields::default())
    }
}

impl Default for Record<5> {
    fn default() -> Self {
        Self(Fields::default())
    }
}

impl Default for Record<6> {
    fn default() -> Self {
        Self(Fields::default())
    }
}

impl crate::feature::Record for Record<3> {
    fn standard_field_count(&self) -> usize {
        3
    }

    fn reference_sequence_name(&self) -> &BStr {
        self.reference_sequence_name()
    }

    fn feature_start(&self) -> io::Result<Position> {
        self.feature_start()
    }

    fn feature_end(&self) -> Option<io::Result<Position>> {
        self.feature_end()
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

    fn other_fields(&self) -> Box<dyn crate::feature::record::OtherFields + '_> {
        Box::new(self.other_fields())
    }
}

impl crate::feature::Record for Record<4> {
    fn standard_field_count(&self) -> usize {
        4
    }

    fn reference_sequence_name(&self) -> &BStr {
        self.reference_sequence_name()
    }

    fn feature_start(&self) -> io::Result<Position> {
        self.feature_start()
    }

    fn feature_end(&self) -> Option<io::Result<Position>> {
        self.feature_end()
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

    fn other_fields(&self) -> Box<dyn crate::feature::record::OtherFields + '_> {
        Box::new(self.other_fields())
    }
}

impl crate::feature::Record for Record<5> {
    fn standard_field_count(&self) -> usize {
        5
    }

    fn reference_sequence_name(&self) -> &BStr {
        self.reference_sequence_name()
    }

    fn feature_start(&self) -> io::Result<Position> {
        self.feature_start()
    }

    fn feature_end(&self) -> Option<io::Result<Position>> {
        self.feature_end()
    }

    fn name(&self) -> Option<Option<&BStr>> {
        Some(self.name())
    }

    fn score(&self) -> Option<io::Result<u16>> {
        Some(self.score())
    }

    fn strand(&self) -> Option<io::Result<Option<Strand>>> {
        None
    }

    fn other_fields(&self) -> Box<dyn crate::feature::record::OtherFields + '_> {
        Box::new(self.other_fields())
    }
}

impl crate::feature::Record for Record<6> {
    fn standard_field_count(&self) -> usize {
        6
    }

    fn reference_sequence_name(&self) -> &BStr {
        self.reference_sequence_name()
    }

    fn feature_start(&self) -> io::Result<Position> {
        self.feature_start()
    }

    fn feature_end(&self) -> Option<io::Result<Position>> {
        self.feature_end()
    }

    fn name(&self) -> Option<Option<&BStr>> {
        Some(self.name())
    }

    fn score(&self) -> Option<io::Result<u16>> {
        Some(self.score())
    }

    fn strand(&self) -> Option<io::Result<Option<Strand>>> {
        Some(self.strand())
    }

    fn other_fields(&self) -> Box<dyn crate::feature::record::OtherFields + '_> {
        Box::new(self.other_fields())
    }
}
