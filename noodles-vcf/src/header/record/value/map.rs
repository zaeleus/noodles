//! VCF header map value.

pub mod alternative_allele;
mod builder;
pub mod contig;
pub mod filter;
pub mod format;
pub mod info;
pub mod meta;
pub mod other;
mod tag;

pub use self::{
    alternative_allele::AlternativeAllele, builder::Builder, contig::Contig, filter::Filter,
    format::Format, info::Info, meta::Meta, other::Other,
};

use std::{
    error,
    fmt::{self, Display},
    str::FromStr,
};

use indexmap::IndexMap;

use crate::header::Number;

type Fields = Vec<(String, String)>;
type OtherFields<S> = IndexMap<tag::Other<S>, String>;

/// An inner VCF header map value.
pub trait Inner: Sized {
    /// The standard tag type.
    type StandardTag: tag::Standard;

    /// The builder type.
    type Builder: builder::Inner<Self>;
}

/// An inner VCF header map value with number and type fields.
pub trait Typed: Inner {
    /// The type type.
    type Type: Display;

    /// Returns the cardinality of the field value.
    fn number(&self) -> Number;

    /// Returns a mutable reference to the number.
    fn number_mut(&mut self) -> &mut Number;

    /// Returns the type of the field value.
    fn ty(&self) -> Self::Type;

    /// Returns a mutable reference to the type.
    fn type_mut(&mut self) -> &mut Self::Type;
}

/// An inner VCF header map value with a description field.
pub trait Described: Inner {
    /// Returns the description.
    fn description(&self) -> &str;

    /// Returns a mutable reference to the description.
    fn description_mut(&mut self) -> &mut String;
}

/// An inner VCF header map value with an IDX field.
pub trait Indexed: Inner {
    /// Returns the index of the ID in the dictionary of strings.
    fn idx(&self) -> Option<usize>;

    /// Returns a mutable reference to the index.
    fn idx_mut(&mut self) -> &mut Option<usize>;
}

/// A VCF header map value.
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
    /// Creates a VCF header map value builder.
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

impl<I> Map<I>
where
    I: Typed,
{
    /// Returns the cardinality of the field value.
    pub fn number(&self) -> Number {
        self.inner.number()
    }

    /// Returns a mutable reference to the number.
    pub fn number_mut(&mut self) -> &mut Number {
        self.inner.number_mut()
    }

    /// Returns the type of the field value.
    pub fn ty(&self) -> I::Type {
        self.inner.ty()
    }

    /// Returns a mutable reference to the type.
    pub fn type_mut(&mut self) -> &mut I::Type {
        self.inner.type_mut()
    }
}

impl<I> Map<I>
where
    I: Described,
{
    /// Returns the description.
    pub fn description(&self) -> &str {
        self.inner.description()
    }

    /// Returns a mutable reference to the description.
    pub fn description_mut(&mut self) -> &mut String {
        self.inner.description_mut()
    }
}

impl<I> Map<I>
where
    I: Indexed,
{
    /// Returns the index of the ID in the dictionary of strings.
    pub fn idx(&self) -> Option<usize> {
        self.inner.idx()
    }

    /// Returns a mutable reference to the index.
    pub fn idx_mut(&mut self) -> &mut Option<usize> {
        self.inner.idx_mut()
    }
}

/// An error returned when a VCF header map value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromFieldsError {
    /// A field is missing.
    MissingField(&'static str),
    /// A tag is duplicated.
    DuplicateTag,
    /// A value is invalid.
    InvalidValue(&'static str),
    /// The actual number does not match the expected number in the reserved definition.
    NumberMismatch(Number, Number),
    /// The actual type does not match the expected type in the reserved definition.
    TypeMismatch(String, String),
}

impl error::Error for TryFromFieldsError {}

impl Display for TryFromFieldsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(tag) => write!(f, "missing field: {tag}"),
            Self::DuplicateTag => "duplicate tag".fmt(f),
            Self::InvalidValue(tag) => write!(f, "invalid value for {tag}"),
            Self::NumberMismatch(actual, expected) => {
                write!(f, "number mismatch: expected {expected}, got {actual}")
            }
            Self::TypeMismatch(actual, expected) => {
                write!(f, "type mismatch: expected {expected}, got {actual}")
            }
        }
    }
}

fn fmt_display_type_fields<T>(f: &mut fmt::Formatter<'_>, number: Number, ty: T) -> fmt::Result
where
    T: Display,
{
    write!(f, ",{}={}", tag::NUMBER, number)?;
    write!(f, ",{}={}", tag::TYPE, ty)?;
    Ok(())
}

fn fmt_display_description_field(f: &mut fmt::Formatter<'_>, description: &str) -> fmt::Result {
    use crate::header::fmt::write_escaped_string;

    write!(f, ",{}=", tag::DESCRIPTION)?;
    write_escaped_string(f, description)?;

    Ok(())
}

fn fmt_display_other_fields<S>(
    f: &mut fmt::Formatter<'_>,
    other_fields: &OtherFields<S>,
) -> fmt::Result {
    use crate::header::fmt::write_escaped_string;

    for (key, value) in other_fields {
        write!(f, ",{key}=")?;
        write_escaped_string(f, value)?;
    }

    Ok(())
}

fn fmt_display_idx_field(f: &mut fmt::Formatter<'_>, idx: usize) -> fmt::Result {
    write!(f, ",{}={}", tag::IDX, idx)
}

fn init_other_fields<S>() -> OtherFields<S> {
    IndexMap::new()
}

fn parse_number(s: &str, number: &mut Option<Number>) -> Result<(), TryFromFieldsError> {
    s.parse()
        .map_err(|_| TryFromFieldsError::InvalidValue(tag::NUMBER))
        .and_then(|value| {
            if number.replace(value).is_none() {
                Ok(())
            } else {
                Err(TryFromFieldsError::DuplicateTag)
            }
        })
}

fn parse_type<T>(s: &str, ty: &mut Option<T>) -> Result<(), TryFromFieldsError>
where
    T: FromStr,
{
    s.parse()
        .map_err(|_| TryFromFieldsError::InvalidValue(tag::TYPE))
        .and_then(|value| {
            if ty.replace(value).is_none() {
                Ok(())
            } else {
                Err(TryFromFieldsError::DuplicateTag)
            }
        })
}

fn parse_description(
    s: String,
    description: &mut Option<String>,
) -> Result<(), TryFromFieldsError> {
    if description.replace(s).is_none() {
        Ok(())
    } else {
        Err(TryFromFieldsError::DuplicateTag)
    }
}

fn parse_idx(s: &str, idx: &mut Option<usize>) -> Result<(), TryFromFieldsError> {
    s.parse()
        .map_err(|_| TryFromFieldsError::InvalidValue(tag::IDX))
        .and_then(|n| {
            if idx.replace(n).is_none() {
                Ok(())
            } else {
                Err(TryFromFieldsError::DuplicateTag)
            }
        })
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
