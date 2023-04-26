//! Inner VCF header FORMAT map value.

pub(crate) mod definition;
mod tag;
mod ty;

pub use self::{tag::Tag, ty::Type};

use core::fmt;

use self::tag::StandardTag;
use super::{
    builder, Described, Fields, Indexed, Inner, Map, OtherFields, TryFromFieldsError, Typed,
};
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

impl TryFrom<Fields> for Map<Format> {
    type Error = TryFromFieldsError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        Self::try_from((FileFormat::default(), fields))
    }
}

impl TryFrom<(FileFormat, Fields)> for Map<Format> {
    type Error = TryFromFieldsError;

    fn try_from((_, fields): (FileFormat, Fields)) -> Result<Self, Self::Error> {
        let mut number = None;
        let mut ty = None;
        let mut description = None;
        let mut idx = None;

        let mut other_fields = OtherFields::new();

        for (key, value) in fields {
            match Tag::from(key) {
                tag::ID => return Err(TryFromFieldsError::DuplicateTag),
                tag::NUMBER => {
                    parse_number(&value).and_then(|v| try_replace(&mut number, tag::NUMBER, v))?
                }
                tag::TYPE => parse_type(&value).and_then(|v| try_replace(&mut ty, tag::TYPE, v))?,
                tag::DESCRIPTION => try_replace(&mut description, tag::DESCRIPTION, value)?,
                tag::IDX => parse_idx(&value).and_then(|v| try_replace(&mut idx, tag::IDX, v))?,
                Tag::Other(t) => try_insert(&mut other_fields, t, value)?,
            }
        }

        let number = number.ok_or(TryFromFieldsError::MissingField("Number"))?;
        let ty = ty.ok_or(TryFromFieldsError::MissingField("Type"))?;
        let description = description.ok_or(TryFromFieldsError::MissingField("Description"))?;

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

fn parse_number(s: &str) -> Result<Number, TryFromFieldsError> {
    s.parse()
        .map_err(|_| TryFromFieldsError::InvalidValue("Number"))
}

fn parse_type(s: &str) -> Result<Type, TryFromFieldsError> {
    s.parse()
        .map_err(|_| TryFromFieldsError::InvalidValue("Type"))
}

fn parse_idx(s: &str) -> Result<usize, TryFromFieldsError> {
    s.parse()
        .map_err(|_| TryFromFieldsError::InvalidValue("IDX"))
}

fn try_replace<T>(option: &mut Option<T>, _: Tag, value: T) -> Result<(), TryFromFieldsError> {
    if option.replace(value).is_none() {
        Ok(())
    } else {
        Err(TryFromFieldsError::DuplicateTag)
    }
}

fn try_insert(
    other_fields: &mut OtherFields<StandardTag>,
    tag: super::tag::Other<StandardTag>,
    value: String,
) -> Result<(), TryFromFieldsError> {
    use indexmap::map::Entry;

    match other_fields.entry(tag) {
        Entry::Vacant(entry) => {
            entry.insert(value);
            Ok(())
        }
        Entry::Occupied(_) => Err(TryFromFieldsError::DuplicateTag),
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
    fn test_try_from_fields_for_map_format() -> Result<(), TryFromFieldsError> {
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
            Err(TryFromFieldsError::MissingField("Number"))
        );

        assert_eq!(
            Map::<Format>::try_from(vec![
                (String::from("Number"), String::from("1")),
                (String::from("Description"), String::from("Genotype")),
            ]),
            Err(TryFromFieldsError::MissingField("Type"))
        );

        assert_eq!(
            Map::<Format>::try_from(vec![
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("String")),
            ]),
            Err(TryFromFieldsError::MissingField("Description"))
        );
    }
}
