use core::fmt;

use indexmap::IndexMap;

use super::{builder, tag, Described, Fields, Indexed, Inner, Map, TryFromFieldsError, Typed};
use crate::header::{
    format::{Key, Type},
    FileFormat, Number,
};

type StandardTag = tag::TypedDescribedIndexed;
type Tag = tag::Tag<StandardTag>;

/// An inner VCF header format map value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Format {
    number: Number,
    ty: Type,
    description: String,
    idx: Option<usize>,
}

impl Inner for Format {
    type Id = Key;
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
    /// use noodles_vcf::header::{
    ///     format::{Key, Type},
    ///     record::value::{map::Format, Map},
    ///     Number,
    /// };
    ///
    /// let map = Map::<Format>::new(Key::Genotype, Number::Count(1), Type::String, "Genotype");
    /// ```
    pub fn new<D>(id: Key, number: Number, ty: Type, description: D) -> Self
    where
        D: Into<String>,
    {
        Self {
            id,
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
        super::fmt_display_prefix(f, self.id())?;
        super::fmt_display_type_fields(f, self.number(), self.ty())?;
        super::fmt_display_description_field(f, self.description())?;
        super::fmt_display_other_fields(f, self.other_fields())?;

        if let Some(idx) = self.idx() {
            super::fmt_display_idx_field(f, idx)?;
        }

        super::fmt_display_suffix(f)?;

        Ok(())
    }
}

impl From<Key> for Map<Format> {
    fn from(key: Key) -> Self {
        use crate::header::format::key;

        let number = key::number(&key).unwrap_or(Number::Count(1));
        let ty = key::ty(&key).unwrap_or(Type::String);
        let description = key::description(&key).map(|s| s.into()).unwrap_or_default();

        Self {
            id: key,
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

    fn try_from((file_format, fields): (FileFormat, Fields)) -> Result<Self, Self::Error> {
        let mut other_fields = super::init_other_fields(fields.len());

        let mut id = None;
        let mut number = None;
        let mut ty = None;
        let mut description = None;
        let mut idx = None;

        for (key, value) in fields {
            match Tag::from(key) {
                Tag::Standard(StandardTag::Id) => super::parse_id(&value, &mut id)?,
                Tag::Standard(StandardTag::Number) => super::parse_number(&value, &mut number)?,
                Tag::Standard(StandardTag::Type) => super::parse_type(&value, &mut ty)?,
                Tag::Standard(StandardTag::Description) => {
                    super::parse_description(value, &mut description)?
                }
                Tag::Standard(StandardTag::Idx) => super::parse_idx(&value, &mut idx)?,
                Tag::Other(t) => super::insert_other_field(&mut other_fields, t, value)?,
            }
        }

        let id = id.ok_or(TryFromFieldsError::MissingField("ID"))?;
        let number = number.ok_or(TryFromFieldsError::MissingField("Number"))?;
        let ty = ty.ok_or(TryFromFieldsError::MissingField("Type"))?;
        let description = description.ok_or(TryFromFieldsError::MissingField("Description"))?;

        if file_format >= FileFormat::new(4, 3) && !matches!(id, Key::Other(_)) {
            validate_type_fields(&id, number, ty)?;
        }

        Ok(Self {
            id,
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

fn validate_type_fields(
    id: &Key,
    actual_number: Number,
    actual_type: Type,
) -> Result<(), TryFromFieldsError> {
    use crate::header::format::key;

    let expected_number = key::number(id).unwrap();

    if actual_number != expected_number {
        return Err(TryFromFieldsError::NumberMismatch);
    }

    let expected_type = key::ty(id).unwrap();

    if actual_type != expected_type {
        return Err(TryFromFieldsError::TypeMismatch);
    }

    Ok(())
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

    #[test]
    fn test_fmt() {
        let map = Map::<Format>::from(Key::Genotype);
        let expected = r#"<ID=GT,Number=1,Type=String,Description="Genotype">"#;
        assert_eq!(map.to_string(), expected);
    }

    #[test]
    fn test_try_from_fields_for_map_format() -> Result<(), TryFromFieldsError> {
        let actual = Map::<Format>::try_from(vec![
            (String::from("ID"), String::from("GT")),
            (String::from("Number"), String::from("1")),
            (String::from("Type"), String::from("String")),
            (String::from("Description"), String::from("Genotype")),
        ])?;

        let expected = Map::<Format>::from(Key::Genotype);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_format_with_missing_fields() {
        assert_eq!(
            Map::<Format>::try_from(vec![
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("String")),
                (String::from("Description"), String::from("Genotype")),
            ]),
            Err(TryFromFieldsError::MissingField("ID"))
        );

        assert_eq!(
            Map::<Format>::try_from(vec![
                (String::from("ID"), String::from("GT")),
                (String::from("Type"), String::from("String")),
                (String::from("Description"), String::from("Genotype")),
            ]),
            Err(TryFromFieldsError::MissingField("Number"))
        );

        assert_eq!(
            Map::<Format>::try_from(vec![
                (String::from("ID"), String::from("GT")),
                (String::from("Number"), String::from("1")),
                (String::from("Description"), String::from("Genotype")),
            ]),
            Err(TryFromFieldsError::MissingField("Type"))
        );

        assert_eq!(
            Map::<Format>::try_from(vec![
                (String::from("ID"), String::from("GT")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("String")),
            ]),
            Err(TryFromFieldsError::MissingField("Description"))
        );
    }
}
