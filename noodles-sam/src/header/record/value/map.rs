//! SAM header record map value.

pub mod builder;
pub mod header;
pub mod program;
pub mod read_group;
pub mod reference_sequence;
pub(crate) mod tag;

pub use self::{
    builder::Builder, header::Header, program::Program, read_group::ReadGroup,
    reference_sequence::ReferenceSequence, tag::Tag,
};

use std::fmt;

use indexmap::IndexMap;

pub(crate) type OtherFields<S> = IndexMap<tag::Other<S>, String>;

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
    pub(crate) inner: I,
    pub(crate) other_fields: OtherFields<I::StandardTag>,
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

fn fmt_display_other_fields<S>(
    f: &mut fmt::Formatter<'_>,
    other_fields: &OtherFields<S>,
) -> fmt::Result {
    const DELIMITER: char = '\t';

    for (key, value) in other_fields {
        write!(f, "{DELIMITER}{key}:{value}")?;
    }

    Ok(())
}
