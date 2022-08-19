//! SAM header record map value builder.

use std::{error, fmt};

use super::{tag, Map, OtherFields};

/// An error returned when a SAM header record map value fails to build.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BuildError {
    /// A required field is missing.
    MissingField(&'static str),
    /// A value is invalid.
    InvalidValue(&'static str),
}

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(tag) => write!(f, "missing field: {}", tag),
            Self::InvalidValue(tag) => write!(f, "invalid value: {}", tag),
        }
    }
}

/// An inner SAM header map value builder.
pub trait Inner<I>: Default
where
    I: super::Inner,
{
    /// Builds the SAM header map value.
    fn build(self) -> Result<I, BuildError>;
}

/// A SAM header record map value builder.
pub struct Builder<I>
where
    I: super::Inner,
{
    pub(crate) inner: I::Builder,
    other_fields: OtherFields<I::StandardTag>,
}

impl<I> Builder<I>
where
    I: super::Inner,
{
    /// Inserts a key-value pair into the other fields.
    pub fn insert(mut self, key: tag::Other<I::StandardTag>, value: String) -> Self {
        self.other_fields.insert(key, value);
        self
    }

    /// Builds a SAM header record map value.
    pub fn build(self) -> Result<Map<I>, BuildError> {
        let inner = self.inner.build()?;

        Ok(Map {
            inner,
            other_fields: self.other_fields,
        })
    }
}

impl<I> Default for Builder<I>
where
    I: super::Inner,
{
    fn default() -> Self {
        Self {
            inner: I::Builder::default(),
            other_fields: OtherFields::new(),
        }
    }
}
