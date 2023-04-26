//! Inner VCF header meta map value.

mod builder;
mod tag;

pub use self::tag::Tag;

use std::fmt::{self, Display};

use self::tag::StandardTag;
use super::{Fields, Inner, Map, OtherFields, TryFromFieldsError};
use crate::header::Number;

/// An inner VCF header meta map value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Meta {
    values: Vec<String>,
}

impl Inner for Meta {
    type StandardTag = StandardTag;
    type Builder = builder::Builder;
}

impl Map<Meta> {
    /// Creates a VCF header meta map value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Meta, Map};
    ///
    /// let id = "Assay";
    /// let map = Map::<Meta>::new(
    ///     vec![String::from("WholeGenome"), String::from("Exome")],
    /// );
    /// ```
    pub fn new(values: Vec<String>) -> Self {
        Self {
            inner: Meta { values },
            other_fields: OtherFields::new(),
        }
    }

    /// Returns the meta values.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Meta, Map};
    /// let values = vec![String::from("WholeGenome"), String::from("Exome")];
    /// let map = Map::<Meta>::new(values.clone());
    /// assert_eq!(map.values(), &values);
    /// ```
    pub fn values(&self) -> &[String] {
        &self.inner.values
    }
}

impl Display for Map<Meta> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        ",Type=String".fmt(f)?;
        write!(f, ",Number={}", Number::Unknown)?;

        ",Values=".fmt(f)?;
        '['.fmt(f)?;

        for (i, value) in self.values().iter().enumerate() {
            if i > 0 {
                ", ".fmt(f)?;
            }

            value.fmt(f)?;
        }

        ']'.fmt(f)?;

        super::fmt_display_other_fields(f, self.other_fields())?;

        Ok(())
    }
}

impl TryFrom<Fields> for Map<Meta> {
    type Error = TryFromFieldsError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut ty = None;
        let mut number = None;
        let mut values = None;

        let mut other_fields = OtherFields::new();

        for (key, value) in fields {
            match Tag::from(key) {
                tag::ID => return Err(TryFromFieldsError::DuplicateTag),
                tag::TYPE => try_replace(&mut ty, tag::TYPE, value)?,
                tag::NUMBER => try_replace(&mut number, tag::NUMBER, value)?,
                tag::VALUES => parse_values(&value, &mut values)?,
                Tag::Other(t) => super::insert_other_field(&mut other_fields, t, value)?,
            }
        }

        let _ = ty.ok_or(TryFromFieldsError::MissingField("Type"))?;
        let _ = number.ok_or(TryFromFieldsError::MissingField("Number"))?;
        let values = values.ok_or(TryFromFieldsError::MissingField("Values"))?;

        Ok(Self {
            inner: Meta { values },
            other_fields,
        })
    }
}

fn parse_values(s: &str, values: &mut Option<Vec<String>>) -> Result<(), TryFromFieldsError> {
    const DELIMITER: char = ',';

    let value = s.split(DELIMITER).map(|t| t.trim().into()).collect();

    if values.replace(value).is_none() {
        Ok(())
    } else {
        Err(TryFromFieldsError::DuplicateTag)
    }
}

fn try_replace<T>(option: &mut Option<T>, _: Tag, value: T) -> Result<(), TryFromFieldsError> {
    if option.replace(value).is_none() {
        Ok(())
    } else {
        Err(TryFromFieldsError::DuplicateTag)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let map = Map::<Meta>::new(vec![String::from("WholeGenome"), String::from("Exome")]);
        let expected = r#",Type=String,Number=.,Values=[WholeGenome, Exome]"#;
        assert_eq!(map.to_string(), expected);
    }

    #[test]
    fn test_try_from_fields_for_map_meta() -> Result<(), TryFromFieldsError> {
        let actual = Map::<Meta>::try_from(vec![
            (String::from("Type"), String::from("String")),
            (String::from("Number"), String::from(".")),
            (String::from("Values"), String::from("WholeGenome, Exome")),
        ])?;

        let expected = Map::<Meta>::new(vec![String::from("WholeGenome"), String::from("Exome")]);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_meta_with_missing_fields() {
        assert_eq!(
            Map::<Meta>::try_from(vec![
                (String::from("Number"), String::from(".")),
                (String::from("Values"), String::from("WholeGenome, Exome")),
            ]),
            Err(TryFromFieldsError::MissingField("Type"))
        );

        assert_eq!(
            Map::<Meta>::try_from(vec![
                (String::from("Type"), String::from("String")),
                (String::from("Values"), String::from("WholeGenome, Exome")),
            ]),
            Err(TryFromFieldsError::MissingField("Number"))
        );

        assert_eq!(
            Map::<Meta>::try_from(vec![
                (String::from("Type"), String::from("String")),
                (String::from("Number"), String::from(".")),
            ]),
            Err(TryFromFieldsError::MissingField("Values"))
        );
    }
}
