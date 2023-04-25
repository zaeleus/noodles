//! Inner VCF header INFO map value.

pub(crate) mod definition;
mod tag;
mod ty;

pub use self::{tag::Tag, ty::Type};

use std::fmt;

use indexmap::IndexMap;

use self::tag::StandardTag;
use super::{builder, Described, Fields, Indexed, Inner, Map, TryFromFieldsError, Typed};
use crate::{
    header::{FileFormat, Number},
    record::info::field::Key,
};

/// An inner VCF header info map value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Info {
    number: Number,
    ty: Type,
    description: String,
    idx: Option<usize>,
}

impl Inner for Info {
    type StandardTag = StandardTag;
    type Builder = builder::TypedDescribedIndexed<Self>;
}

impl Typed for Info {
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

impl Described for Info {
    fn description(&self) -> &str {
        &self.description
    }

    fn description_mut(&mut self) -> &mut String {
        &mut self.description
    }
}

impl Indexed for Info {
    fn idx(&self) -> Option<usize> {
        self.idx
    }

    fn idx_mut(&mut self) -> &mut Option<usize> {
        &mut self.idx
    }
}

impl Map<Info> {
    /// Creates a VCF header info map value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::{record::value::{map::{info::Type, Info}, Map}, Number},
    ///     record::info::field::key,
    /// };
    ///
    /// let id = key::SAMPLES_WITH_DATA_COUNT;
    /// let map = Map::<Info>::new(
    ///     Number::Count(1),
    ///     Type::Integer,
    ///     "Number of samples with data",
    /// );
    /// ```
    pub fn new<D>(number: Number, ty: Type, description: D) -> Self
    where
        D: Into<String>,
    {
        Self {
            inner: Info {
                number,
                ty,
                description: description.into(),
                idx: None,
            },
            other_fields: IndexMap::new(),
        }
    }
}

impl fmt::Display for Map<Info> {
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

impl From<&Key> for Map<Info> {
    fn from(key: &Key) -> Self {
        Self::from((FileFormat::default(), key))
    }
}

impl From<(FileFormat, &Key)> for Map<Info> {
    fn from((file_format, key): (FileFormat, &Key)) -> Self {
        let (number, ty, description) =
            definition::definition(file_format, key).unwrap_or_default();

        Self {
            inner: Info {
                number,
                ty,
                description: description.into(),
                idx: None,
            },
            other_fields: IndexMap::new(),
        }
    }
}

impl TryFrom<Fields> for Map<Info> {
    type Error = TryFromFieldsError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        Self::try_from((FileFormat::default(), fields))
    }
}

impl TryFrom<(FileFormat, Fields)> for Map<Info> {
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
            inner: Info {
                number,
                ty,
                description,
                idx,
            },
            other_fields,
        })
    }
}

impl builder::Inner<Info> for builder::TypedDescribedIndexed<Info> {
    fn build(self) -> Result<Info, builder::BuildError> {
        let number = self
            .number
            .ok_or(builder::BuildError::MissingField("Number"))?;

        let ty = self.ty.ok_or(builder::BuildError::MissingField("Type"))?;

        let description = self
            .description
            .ok_or(builder::BuildError::MissingField("Description"))?;

        Ok(Info {
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
    use crate::record::info::field::key;

    #[test]
    fn test_fmt() {
        let map = Map::<Info>::from(&key::SAMPLES_WITH_DATA_COUNT);
        let expected = r#",Number=1,Type=Integer,Description="Number of samples with data""#;
        assert_eq!(map.to_string(), expected);
    }

    #[test]
    fn test_try_from_fields_for_map_info() -> Result<(), TryFromFieldsError> {
        let actual = Map::<Info>::try_from(vec![
            (String::from("Number"), String::from("1")),
            (String::from("Type"), String::from("Integer")),
            (
                String::from("Description"),
                String::from("Number of samples with data"),
            ),
        ])?;

        let expected = Map::<Info>::from(&key::SAMPLES_WITH_DATA_COUNT);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_info_with_missing_fields() {
        assert_eq!(
            Map::<Info>::try_from(vec![
                (String::from("Type"), String::from("Integer")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data")
                ),
            ]),
            Err(TryFromFieldsError::MissingField("Number"))
        );

        assert_eq!(
            Map::<Info>::try_from(vec![
                (String::from("Number"), String::from("1")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data")
                ),
            ]),
            Err(TryFromFieldsError::MissingField("Type"))
        );

        assert_eq!(
            Map::<Info>::try_from(vec![
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("Integer")),
            ]),
            Err(TryFromFieldsError::MissingField("Description"))
        );
    }
}
