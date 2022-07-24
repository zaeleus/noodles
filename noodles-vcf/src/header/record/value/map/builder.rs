use std::{error, fmt};

use super::{Map, OtherFields};
use crate::header::Number;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BuildError {
    MissingField(&'static str),
}

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(tag) => write!(f, "missing field: {}", tag),
        }
    }
}

pub trait Inner<I>: Default
where
    I: super::Inner,
{
    fn build(self) -> Result<I, BuildError>;
}

pub trait Typed<I>: Inner<I>
where
    I: super::Typed,
{
    fn set_number(self, number: Number) -> Self;
    fn set_type(self, ty: I::Type) -> Self;
}

pub trait Described<I>: Inner<I>
where
    I: super::Described,
{
    fn set_description<D>(self, description: D) -> Self
    where
        D: Into<String>;
}

pub trait Indexed<I>: Inner<I>
where
    I: super::Indexed,
{
    fn set_idx(self, idx: usize) -> Self;
}

/// A VCF header record map value builder.
pub struct Builder<I>
where
    I: super::Inner,
{
    id: Option<I::Id>,
    pub(super) inner: I::Builder,
    other_fields: OtherFields<I::StandardTag>,
}

impl<I> Builder<I>
where
    I: super::Inner,
{
    /// Sets the ID.
    pub fn set_id(mut self, id: I::Id) -> Self {
        self.id = Some(id);
        self
    }

    /// Builds a VCF header record map value.
    pub fn build(self) -> Result<Map<I>, BuildError> {
        let id = self.id.ok_or(BuildError::MissingField("ID"))?;
        let inner = self.inner.build()?;

        Ok(Map {
            id,
            inner,
            other_fields: self.other_fields,
        })
    }
}

impl<I> Builder<I>
where
    I: super::Typed,
    I::Builder: Typed<I>,
{
    /// Sets the number.
    pub fn set_number(mut self, number: Number) -> Self {
        self.inner = self.inner.set_number(number);
        self
    }

    /// Sets the type.
    pub fn set_type(mut self, ty: I::Type) -> Self {
        self.inner = self.inner.set_type(ty);
        self
    }
}

impl<I> Builder<I>
where
    I: super::Described,
    I::Builder: Described<I>,
{
    /// Sets the description.
    pub fn set_description<D>(mut self, description: D) -> Self
    where
        D: Into<String>,
    {
        self.inner = self.inner.set_description(description);
        self
    }
}

impl<I> Builder<I>
where
    I: super::Indexed,
    I::Builder: Indexed<I>,
{
    /// Sets the index.
    pub fn set_idx(mut self, idx: usize) -> Self {
        self.inner = self.inner.set_idx(idx);
        self
    }
}

impl<I> Default for Builder<I>
where
    I: super::Inner,
{
    fn default() -> Self {
        Self {
            id: None,
            inner: I::Builder::default(),
            other_fields: OtherFields::new(),
        }
    }
}

#[derive(Default)]
pub struct Identity;

impl<I> Inner<I> for Identity
where
    I: Default + super::Inner,
{
    fn build(self) -> Result<I, BuildError> {
        Ok(I::default())
    }
}

pub struct TypedDescribedIndexed<I>
where
    I: super::Typed + super::Described + super::Indexed,
{
    pub(super) number: Option<Number>,
    pub(super) ty: Option<I::Type>,
    pub(super) description: Option<String>,
    pub(super) idx: Option<usize>,
}

impl<I> Default for TypedDescribedIndexed<I>
where
    I: super::Typed + super::Described + super::Indexed,
{
    fn default() -> Self {
        Self {
            number: None,
            ty: None,
            description: None,
            idx: None,
        }
    }
}

#[derive(Default)]
pub struct DescribedIndexed {
    pub(super) description: Option<String>,
    pub(super) idx: Option<usize>,
}
