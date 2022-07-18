use std::fmt;

use indexmap::IndexMap;

use super::{Described, Fields, Indexed, Inner, Map, TryFromFieldsError, Typed};
use crate::header::{
    info::{Key, Type},
    FileFormat, Number,
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
    type Id = Key;
}

impl Typed for Info {
    type Type = Type;

    fn number(&self) -> Number {
        self.number
    }

    fn ty(&self) -> Self::Type {
        self.ty
    }
}

impl Described for Info {
    fn description(&self) -> &str {
        &self.description
    }
}

impl Indexed for Info {
    fn idx(&self) -> Option<usize> {
        self.idx
    }
}

impl Map<Info> {
    /// Creates a VCF header info map value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{
    ///     info::{Key, Type},
    ///     record::value::{map::Info, Map},
    ///     Number,
    /// };
    ///
    /// let map = Map::<Info>::new(
    ///     Key::SamplesWithDataCount,
    ///     Number::Count(1),
    ///     Type::Integer,
    ///     "Number of samples with data",
    /// );
    /// ```
    pub fn new<D>(id: Key, number: Number, ty: Type, description: D) -> Self
    where
        D: Into<String>,
    {
        Self {
            id,
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

impl From<Key> for Map<Info> {
    fn from(key: Key) -> Self {
        use crate::header::info::key;

        let number = key::number(&key).unwrap_or(Number::Count(1));
        let ty = key::ty(&key).unwrap_or(Type::String);
        let description = key::description(&key).map(|s| s.into()).unwrap_or_default();

        Self {
            id: key,
            inner: Info {
                number,
                ty,
                description,
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

    fn try_from((file_format, fields): (FileFormat, Fields)) -> Result<Self, Self::Error> {
        let mut other_fields = super::init_other_fields(fields.len());

        let mut id = None;
        let mut number = None;
        let mut ty = None;
        let mut description = None;
        let mut idx = None;

        for (key, value) in fields {
            match key.as_str() {
                "ID" => super::parse_id(&value, &mut id)?,
                "Number" => super::parse_number(&value, &mut number)?,
                "Type" => super::parse_type(&value, &mut ty)?,
                "Description" => super::parse_description(value, &mut description)?,
                "IDX" => super::parse_idx(&value, &mut idx)?,
                _ => super::insert_other_field(&mut other_fields, key, value)?,
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

fn validate_type_fields(
    id: &Key,
    actual_number: Number,
    actual_type: Type,
) -> Result<(), TryFromFieldsError> {
    use crate::header::info::key;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let map = Map::<Info>::from(Key::SamplesWithDataCount);
        let expected = r#"<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">"#;
        assert_eq!(map.to_string(), expected);
    }

    #[test]
    fn test_try_from_fields_for_map_info() -> Result<(), TryFromFieldsError> {
        let actual = Map::<Info>::try_from(vec![
            (String::from("ID"), String::from("NS")),
            (String::from("Number"), String::from("1")),
            (String::from("Type"), String::from("Integer")),
            (
                String::from("Description"),
                String::from("Number of samples with data"),
            ),
        ])?;

        let expected = Map::<Info>::from(Key::SamplesWithDataCount);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_info_with_missing_fields() {
        assert_eq!(
            Map::<Info>::try_from(vec![
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("Integer")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data")
                ),
            ]),
            Err(TryFromFieldsError::MissingField("ID"))
        );

        assert_eq!(
            Map::<Info>::try_from(vec![
                (String::from("ID"), String::from("NS")),
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
                (String::from("ID"), String::from("NS")),
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
                (String::from("ID"), String::from("NS")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("Integer")),
            ]),
            Err(TryFromFieldsError::MissingField("Description"))
        );
    }
}
