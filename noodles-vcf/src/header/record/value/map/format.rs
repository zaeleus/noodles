//! Inner VCF header FORMAT map value.

pub(crate) mod definition;
mod ty;

pub use self::ty::Type;

use core::fmt;

use indexmap::IndexMap;

use super::{builder, tag, Described, Fields, Indexed, Inner, Map, TryFromFieldsError, Typed};
use crate::{
    header::{FileFormat, Number},
    record::genotypes::keys::Key,
};

type StandardTag = tag::TypedDescribedIndexed;

/// A VCF header format map tag.
pub type Tag = tag::Tag<StandardTag>;

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
            other_fields: IndexMap::new(),
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
        let number = definition::number(key).unwrap_or_default();
        let ty = definition::ty(key).unwrap_or(Type::String);
        let description = definition::description(key)
            .map(|s| s.into())
            .unwrap_or_default();

        Self {
            inner: Format {
                number,
                ty,
                description,
                idx: None,
            },
            other_fields: IndexMap::new(),
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
        let mut other_fields = super::init_other_fields();

        let mut number = None;
        let mut ty = None;
        let mut description = None;
        let mut idx = None;

        for (key, value) in fields {
            match Tag::from(key) {
                Tag::Standard(StandardTag::Id) => return Err(TryFromFieldsError::DuplicateTag),
                Tag::Standard(StandardTag::Number) => super::parse_number(&value, &mut number)?,
                Tag::Standard(StandardTag::Type) => super::parse_type(&value, &mut ty)?,
                Tag::Standard(StandardTag::Description) => {
                    super::parse_description(value, &mut description)?
                }
                Tag::Standard(StandardTag::Idx) => super::parse_idx(&value, &mut idx)?,
                Tag::Other(t) => super::insert_other_field(&mut other_fields, t, value)?,
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
