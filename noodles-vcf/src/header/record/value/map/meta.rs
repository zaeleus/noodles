//! Inner VCF header meta map value.

mod builder;
pub(crate) mod tag;

pub use self::tag::Tag;

use std::{error, fmt};

use self::tag::StandardTag;
use super::{Fields, Inner, Map, OtherFields};
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

impl fmt::Display for Map<Meta> {
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

/// An error returned when a raw META record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is missing.
    MissingField(Tag),
    /// A tag is duplicated.
    DuplicateTag(Tag),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(tag) => write!(f, "missing field: {tag}"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

impl TryFrom<Fields> for Map<Meta> {
    type Error = ParseError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut ty = None;
        let mut number = None;
        let mut values = None;

        let mut other_fields = OtherFields::new();

        for (key, value) in fields {
            match Tag::from(key) {
                tag::ID => return Err(ParseError::DuplicateTag(tag::ID)),
                tag::TYPE => try_replace(&mut ty, tag::TYPE, value)?,
                tag::NUMBER => try_replace(&mut number, tag::NUMBER, value)?,
                tag::VALUES => {
                    let v = parse_values(&value);
                    try_replace(&mut values, tag::VALUES, v)?;
                }
                Tag::Other(t) => try_insert(&mut other_fields, t, value)?,
            }
        }

        let _ = ty.ok_or(ParseError::MissingField(tag::TYPE))?;
        let _ = number.ok_or(ParseError::MissingField(tag::NUMBER))?;
        let values = values.ok_or(ParseError::MissingField(tag::VALUES))?;

        Ok(Self {
            inner: Meta { values },
            other_fields,
        })
    }
}

fn parse_values(s: &str) -> Vec<String> {
    const DELIMITER: char = ',';
    s.split(DELIMITER).map(|t| t.trim().into()).collect()
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
    fn test_try_from_fields_for_map_meta() -> Result<(), ParseError> {
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
            Err(ParseError::MissingField(tag::TYPE))
        );

        assert_eq!(
            Map::<Meta>::try_from(vec![
                (String::from("Type"), String::from("String")),
                (String::from("Values"), String::from("WholeGenome, Exome")),
            ]),
            Err(ParseError::MissingField(tag::NUMBER))
        );

        assert_eq!(
            Map::<Meta>::try_from(vec![
                (String::from("Type"), String::from("String")),
                (String::from("Number"), String::from(".")),
            ]),
            Err(ParseError::MissingField(tag::VALUES))
        );
    }
}
