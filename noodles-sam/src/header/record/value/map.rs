//! SAM header record map value.

pub mod builder;
pub mod header;
mod program;
pub mod read_group;
pub mod reference_sequence;
mod tag;

pub use self::{
    builder::Builder, header::Header, program::Program, read_group::ReadGroup,
    reference_sequence::ReferenceSequence, tag::Tag,
};

use std::{error, fmt};

use indexmap::IndexMap;

type Fields = Vec<(String, String)>;
type OtherFields<S> = IndexMap<tag::Other<S>, String>;

/// An inner SAM header record map value.
pub trait Inner: Sized {
    /// The standard tag type.
    type StandardTag: tag::Standard;

    /// The builder type.
    type Builder: builder::Inner<Self>;
}

/// A SAM header record map value.
// TODO: I'm not actually sure what's causing this lint warning.
#[allow(clippy::derive_partial_eq_without_eq)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Map<I>
where
    I: Inner,
{
    inner: I,
    other_fields: OtherFields<I::StandardTag>,
}

impl<I> Map<I>
where
    I: Inner,
{
    /// Creates a SAM header record map value.
    pub fn builder() -> Builder<I> {
        Builder::default()
    }

    /// Returns the nonstandard fields in the map.
    pub fn other_fields(&self) -> &OtherFields<I::StandardTag> {
        &self.other_fields
    }
}

impl<I> Default for Map<I>
where
    I: Inner + Default,
{
    fn default() -> Self {
        Self {
            inner: I::default(),
            other_fields: OtherFields::new(),
        }
    }
}

/// An error returned when a SAM header record map value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromFieldsError {
    /// A field is missing.
    MissingField(&'static str),
    /// A tag is invalid.
    InvalidTag,
    /// A tag is duplicated.
    DuplicateTag,
    /// A value is invalid.
    InvalidValue(&'static str),
}

impl error::Error for TryFromFieldsError {}

impl fmt::Display for TryFromFieldsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(tag) => write!(f, "missing field: {}", tag),
            Self::InvalidTag => "invalid tag".fmt(f),
            Self::DuplicateTag => "duplicate tag".fmt(f),
            Self::InvalidValue(tag) => write!(f, "invalid value for {}", tag),
        }
    }
}

fn fmt_display_other_fields<S>(
    f: &mut fmt::Formatter<'_>,
    other_fields: &OtherFields<S>,
) -> fmt::Result {
    for (key, value) in other_fields {
        write!(f, "\t{}:{}", key, value)?;
    }

    Ok(())
}

fn init_other_fields<S>(len: usize) -> OtherFields<S> {
    let len = len.checked_sub(1).unwrap_or_default();
    IndexMap::with_capacity(len)
}

fn insert_other_field<S>(
    other_fields: &mut OtherFields<S>,
    key: tag::Other<S>,
    value: String,
) -> Result<(), TryFromFieldsError> {
    if other_fields.insert(key, value).is_none() {
        Ok(())
    } else {
        Err(TryFromFieldsError::DuplicateTag)
    }
}
