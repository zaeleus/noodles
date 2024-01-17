//! SAM header record map value builder.

use std::{error, fmt};

use bstr::BString;

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
            Self::MissingField(tag) => write!(f, "missing field: {tag}"),
            Self::InvalidValue(tag) => write!(f, "invalid value: {tag}"),
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
    ///
    /// # Examples
    ///
    /// ```
    /// use bstr::BString;
    /// use noodles_sam::header::record::value::{map::{Header, Tag}, Map};
    ///
    /// let nd = match Tag::try_from([b'n', b'd']) {
    ///     Ok(Tag::Other(tag)) => tag,
    ///     _ => unreachable!(),
    /// };
    ///
    /// let header = Map::<Header>::builder().insert(nd, "noodles").build()?;
    ///
    /// assert_eq!(header.other_fields().get(b"nd"), Some(&BString::from("noodles")));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn insert<V>(mut self, key: tag::Other<I::StandardTag>, value: V) -> Self
    where
        V: Into<BString>,
    {
        self.other_fields.insert(key, value.into());
        self
    }

    /// Builds a SAM header record map value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::Header, Map};
    /// let header = Map::<Header>::builder().build()?;
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
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
