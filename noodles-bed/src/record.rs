pub(crate) mod fields;
mod other_fields;

use std::{fmt, io};

use bstr::ByteSlice;
use noodles_core::Position;

use self::fields::Fields;
pub use self::other_fields::OtherFields;
use crate::feature::record_buf::Strand;

/// A BED record.
#[derive(Clone, Eq, PartialEq)]
pub struct Record<const N: usize>(pub(crate) Fields<N>);

impl<const N: usize> Record<N> {
    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &[u8] {
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
    pub fn name(&self) -> Option<&[u8]> {
        self.0.name()
    }

    /// Returns the score.
    pub fn score(&self) -> Option<io::Result<u16>> {
        self.0.score()
    }

    /// Returns the strand.
    pub fn strand(&self) -> Option<io::Result<Option<Strand>>> {
        self.0.strand()
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

    /// Returns the number of standard fields.
    pub const fn standard_field_count(&self) -> usize {
        N
    }
}

impl<const N: usize> fmt::Debug for Record<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Record")
            .field(
                "reference_sequence_name",
                &self.reference_sequence_name().as_bstr(),
            )
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

impl<const N: usize> crate::feature::Record for Record<N> {
    fn reference_sequence_name(&self) -> &[u8] {
        self.reference_sequence_name()
    }

    fn feature_start(&self) -> io::Result<Position> {
        self.feature_start()
    }

    fn feature_end(&self) -> Option<io::Result<Position>> {
        self.feature_end()
    }
}
