//! Inner VCF header FORMAT map value.

pub(crate) mod definition;
pub(crate) mod tag;
mod ty;

pub use self::{tag::Tag, ty::Type};

use std::{error, fmt, num};

use self::tag::StandardTag;
use super::{builder, Described, Fields, Indexed, Inner, Map, OtherFields, Typed};
use crate::{
    header::{FileFormat, Number},
    record::genotypes::keys::Key,
};

/// An inner VCF header format map value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Format {
    number: Number,
    ty: Type,
    description: String,
    idx: Option<usize>,
}

impl Inner for Format {
    type StandardTag = StandardTag;
    type Builder = builder::TypedDescribedIndexed<Self>;
}

impl Typed for Format {
    type Type = Type;

    fn number(&self) -> Number {
        self.number
    }

    fn number_mut(&mut self) -> &mut Number {
        &mut self.number
    }

    fn ty(&self) -> Self::Type {
        self.ty
    }

    fn type_mut(&mut self) -> &mut Self::Type {
        &mut self.ty
    }
}

impl Described for Format {
    fn description(&self) -> &str {
        &self.description
    }

    fn description_mut(&mut self) -> &mut String {
        &mut self.description
    }
}

impl Indexed for Format {
    fn idx(&self) -> Option<usize> {
        self.idx
    }

    fn idx_mut(&mut self) -> &mut Option<usize> {
        &mut self.idx
    }
}

impl Map<Format> {
    /// Creates a VCF header format map value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::{record::value::{map::{format::Type, Format}, Map}, Number},
    ///     record::genotypes::keys::key,
    /// };
    ///
    /// let id = key::GENOTYPE;
    /// let map = Map::<Format>::new(Number::Count(1), Type::String, "Genotype");
    /// ```
    pub fn new<D>(number: Number, ty: Type, description: D) -> Self
    where
        D: Into<String>,
    {
        Self {
            inner: Format {
                number,
                ty,
                description: description.into(),
                idx: None,
            },
            other_fields: OtherFields::new(),
        }
    }
}

impl fmt::Display for Map<Format> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        super::fmt_display_type_fields(f, self.number(), self.ty())?;
        super::fmt_display_description_field(f, self.description())?;
        super::fmt_display_other_fields(f, self.other_fields())?;

        if let Some(idx) = self.idx() {
            super::fmt_display_idx_field(f, idx)?;
        }

        Ok(())
    }
}

impl From<&Key> for Map<Format> {
    fn from(key: &Key) -> Self {
        Self::from((FileFormat::default(), key))
    }
}

impl From<(FileFormat, &Key)> for Map<Format> {
    fn from((file_format, key): (FileFormat, &Key)) -> Self {
        let (number, ty, description) =
            definition::definition(file_format, key).unwrap_or_default();

        Self {
            inner: Format {
                number,
                ty,
                description: description.into(),
                idx: None,
            },
            other_fields: OtherFields::new(),
        }
    }
}

/// An error returned when a raw FORMAT record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is missing.
    MissingField(Tag),
    /// A tag is duplicated.
    DuplicateTag(Tag),
    /// The ID is invalid.
    InvalidId(crate::record::genotypes::keys::key::ParseError),
    /// The number is invalid.
    InvalidNumber(crate::header::number::ParseError),
    /// The type is invalid.
    InvalidType(ty::ParseError),
    /// The IDX is invalid.
    InvalidIdx(num::ParseIntError),
    /// The number for the given ID does not match its reserved type definition.
    NumberMismatch {
        /// The actual number.
        actual: Number,
        /// The expected number.
        expected: Number,
    },
    /// The type for the given ID does not match its reserved type definition.
    TypeMismatch {
        /// The actual type.
        actual: Type,
        /// The expected type.
        expected: Type,
    },
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidId(e) => Some(e),
            Self::InvalidNumber(e) => Some(e),
            Self::InvalidType(e) => Some(e),
            Self::InvalidIdx(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(tag) => write!(f, "missing field: {tag}"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
            Self::InvalidId(_) => write!(f, "invalid ID"),
            Self::InvalidNumber(_) => write!(f, "invalid number"),
            Self::InvalidType(_) => write!(f, "invalid type"),
            Self::InvalidIdx(_) => write!(f, "invalid IDX"),
            Self::NumberMismatch { actual, expected } => {
                write!(f, "number mismatch: expected {expected}, got {actual}")
            }
            Self::TypeMismatch { actual, expected } => {
                write!(f, "type mismatch: expected {expected}, got {actual}")
            }
        }
    }
}

impl TryFrom<Fields> for Map<Format> {
    type Error = ParseError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        Self::try_from((FileFormat::default(), fields))
    }
}

impl TryFrom<(FileFormat, Fields)> for Map<Format> {
    type Error = ParseError;

    fn try_from((_, fields): (FileFormat, Fields)) -> Result<Self, Self::Error> {
        let mut number = None;
        let mut ty = None;
        let mut description = None;
        let mut idx = None;

        let mut other_fields = OtherFields::new();

        for (key, value) in fields {
            match Tag::from(key) {
                tag::ID => return Err(ParseError::DuplicateTag(tag::ID)),
                tag::NUMBER => {
                    parse_number(&value).and_then(|v| try_replace(&mut number, tag::NUMBER, v))?
                }
                tag::TYPE => parse_type(&value).and_then(|v| try_replace(&mut ty, tag::TYPE, v))?,
                tag::DESCRIPTION => try_replace(&mut description, tag::DESCRIPTION, value)?,
                tag::IDX => parse_idx(&value).and_then(|v| try_replace(&mut idx, tag::IDX, v))?,
                Tag::Other(t) => try_insert(&mut other_fields, t, value)?,
            }
        }

        let number = number.ok_or(ParseError::MissingField(tag::NUMBER))?;
        let ty = ty.ok_or(ParseError::MissingField(tag::TYPE))?;
        let description = description.ok_or(ParseError::MissingField(tag::DESCRIPTION))?;

        Ok(Self {
            inner: Format {
                number,
                ty,
                description,
                idx,
            },
            other_fields,
        })
    }
}

fn parse_number(s: &str) -> Result<Number, ParseError> {
    s.parse().map_err(ParseError::InvalidNumber)
}

fn parse_type(s: &str) -> Result<Type, ParseError> {
    s.parse().map_err(ParseError::InvalidType)
}

fn parse_idx(s: &str) -> Result<usize, ParseError> {
    s.parse().map_err(ParseError::InvalidIdx)
}

fn try_replace<T>(option: &mut Option<T>, tag: Tag, value: T) -> Result<(), ParseError> {
    if option.replace(value).is_none() {
        Ok(())
    } else {
        Err(ParseError::DuplicateTag(tag))
    }
}

fn try_insert(
    other_fields: &mut OtherFields<StandardTag>,
    tag: super::tag::Other<StandardTag>,
    value: String,
) -> Result<(), ParseError> {
    use indexmap::map::Entry;

    match other_fields.entry(tag) {
        Entry::Vacant(entry) => {
            entry.insert(value);
            Ok(())
        }
        Entry::Occupied(entry) => {
            let (t, _) = entry.remove_entry();
            Err(ParseError::DuplicateTag(Tag::Other(t)))
        }
    }
}

impl builder::Inner<Format> for builder::TypedDescribedIndexed<Format> {
    fn build(self) -> Result<Format, builder::BuildError> {
        let number = self
            .number
            .ok_or(builder::BuildError::MissingField("Number"))?;

        let ty = self.ty.ok_or(builder::BuildError::MissingField("Type"))?;

        let description = self
            .description
            .ok_or(builder::BuildError::MissingField("Description"))?;

        Ok(Format {
            number,
            ty,
            description,
            idx: self.idx,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::genotypes::keys::key;

    #[test]
    fn test_fmt() {
        let map = Map::<Format>::from(&key::GENOTYPE);
        let expected = r#",Number=1,Type=String,Description="Genotype""#;
        assert_eq!(map.to_string(), expected);
    }

    #[test]
    fn test_try_from_fields_for_map_format() -> Result<(), ParseError> {
        let actual = Map::<Format>::try_from(vec![
            (String::from("Number"), String::from("1")),
            (String::from("Type"), String::from("String")),
            (String::from("Description"), String::from("Genotype")),
        ])?;

        let expected = Map::<Format>::from(&key::GENOTYPE);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_format_with_missing_fields() {
        assert_eq!(
            Map::<Format>::try_from(vec![
                (String::from("Type"), String::from("String")),
                (String::from("Description"), String::from("Genotype")),
            ]),
            Err(ParseError::MissingField(tag::NUMBER))
        );

        assert_eq!(
            Map::<Format>::try_from(vec![
                (String::from("Number"), String::from("1")),
                (String::from("Description"), String::from("Genotype")),
            ]),
            Err(ParseError::MissingField(tag::TYPE))
        );

        assert_eq!(
            Map::<Format>::try_from(vec![
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("String")),
            ]),
            Err(ParseError::MissingField(tag::DESCRIPTION))
        );
    }
}
