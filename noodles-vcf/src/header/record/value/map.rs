//! VCF header map value.

mod alternative_allele;
mod contig;
mod filter;
mod format;
mod info;
mod meta;
mod other;

pub use self::{
    alternative_allele::AlternativeAllele, contig::Contig, filter::Filter, format::Format,
    info::Info, meta::Meta, other::Other,
};

use std::{
    error,
    fmt::{self, Display},
    str::FromStr,
};

use indexmap::IndexMap;

use crate::header::Number;

type Fields = Vec<(String, String)>;

/// An inner VCF header map value.
pub trait Inner {
    /// The ID type.
    type Id: Display;
}

/// An inner VCF header map value with number and type fields.
pub trait Typed: Inner {
    /// The type type.
    type Type: Display;

    /// Returns the cardinality of the field value.
    fn number(&self) -> Number;

    /// Returns the type of the field value.
    fn ty(&self) -> Self::Type;
}

/// An inner VCF header map value with a description field.
pub trait Described: Inner {
    /// Returns the description.
    fn description(&self) -> &str;
}

/// An inner VCF header map value with an IDX field.
pub trait Indexed: Inner {
    /// Returns the index of the ID in the dictionary of strings.
    fn idx(&self) -> Option<usize>;
}

/// A VCF header map value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Map<I>
where
    I: Inner,
{
    id: I::Id,
    inner: I,
    other_fields: IndexMap<String, String>,
}

impl<I> Map<I>
where
    I: Inner,
{
    /// Returns the ID.
    pub fn id(&self) -> &I::Id {
        &self.id
    }

    /// Returns the nonstandard fields in the map.
    pub fn other_fields(&self) -> &IndexMap<String, String> {
        &self.other_fields
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

    /// Returns the type of the field value.
    pub fn ty(&self) -> I::Type {
        self.inner.ty()
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
}

impl<I> Map<I>
where
    I: Indexed,
{
    /// Returns the index of the ID in the dictionary of strings.
    pub fn idx(&self) -> Option<usize> {
        self.inner.idx()
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
    NumberMismatch,
    /// The actual type does not match the expected type in the reserved definition.
    TypeMismatch,
}

impl error::Error for TryFromFieldsError {}

impl Display for TryFromFieldsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(tag) => write!(f, "missing field: {}", tag),
            Self::DuplicateTag => "duplicate tag".fmt(f),
            Self::InvalidValue(tag) => write!(f, "invalid value for {}", tag),
            Self::NumberMismatch => "number mismatch".fmt(f),
            Self::TypeMismatch => "type mismatch".fmt(f),
        }
    }
}

fn fmt_display_prefix<I>(f: &mut fmt::Formatter<'_>, id: I) -> fmt::Result
where
    I: Display,
{
    write!(f, "<ID={}", id)
}

fn fmt_display_type_fields<T>(f: &mut fmt::Formatter<'_>, number: Number, ty: T) -> fmt::Result
where
    T: Display,
{
    write!(f, ",Number={}", number)?;
    write!(f, ",Type={}", ty)?;
    Ok(())
}

fn fmt_display_description_field(f: &mut fmt::Formatter<'_>, description: &str) -> fmt::Result {
    use crate::header::fmt::write_escaped_string;

    ",Description=".fmt(f)?;
    write_escaped_string(f, description)?;

    Ok(())
}

fn fmt_display_other_fields(
    f: &mut fmt::Formatter<'_>,
    other_fields: &IndexMap<String, String>,
) -> fmt::Result {
    use crate::header::fmt::write_escaped_string;

    for (key, value) in other_fields {
        write!(f, ",{}=", key)?;
        write_escaped_string(f, value)?;
    }

    Ok(())
}

fn fmt_display_idx_field(f: &mut fmt::Formatter<'_>, idx: usize) -> fmt::Result {
    write!(f, ",IDX={}", idx)
}

fn fmt_display_suffix(f: &mut fmt::Formatter<'_>) -> fmt::Result {
    '>'.fmt(f)
}

fn init_other_fields(len: usize) -> IndexMap<String, String> {
    let len = len.checked_sub(1).unwrap_or_default();
    IndexMap::with_capacity(len)
}

fn parse_id<I>(s: &str, id: &mut Option<I>) -> Result<(), TryFromFieldsError>
where
    I: FromStr,
{
    s.parse()
        .map_err(|_| TryFromFieldsError::InvalidValue("ID"))
        .and_then(|value| {
            if id.replace(value).is_none() {
                Ok(())
            } else {
                Err(TryFromFieldsError::DuplicateTag)
            }
        })
}

fn parse_number(s: &str, number: &mut Option<Number>) -> Result<(), TryFromFieldsError> {
    s.parse()
        .map_err(|_| TryFromFieldsError::InvalidValue("Number"))
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
        .map_err(|_| TryFromFieldsError::InvalidValue("Type"))
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
        .map_err(|_| TryFromFieldsError::InvalidValue("IDX"))
        .and_then(|n| {
            if idx.replace(n).is_none() {
                Ok(())
            } else {
                Err(TryFromFieldsError::DuplicateTag)
            }
        })
}

fn insert_other_field(
    other_fields: &mut IndexMap<String, String>,
    key: String,
    value: String,
) -> Result<(), TryFromFieldsError> {
    if other_fields.insert(key, value).is_none() {
        Ok(())
    } else {
        Err(TryFromFieldsError::DuplicateTag)
    }
}
