//! GFF record attributes.

pub mod field;

use std::{
    error, fmt,
    ops::{Deref, DerefMut},
    str::{self, FromStr},
};

use indexmap::IndexMap;

use self::field::{Tag, Value};

const DELIMITER: char = ';';

/// GFF record attributes.
///
/// Attributes are extra data attached to a GFF record. They are represented as a typed map, where
/// each key ([`Tag`]) is associated with a typed [`Value`].
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Attributes(IndexMap<Tag, Value>);

impl Deref for Attributes {
    type Target = IndexMap<Tag, Value>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Attributes {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl fmt::Display for Attributes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use self::field::field_fmt;

        for (i, field) in self.iter().enumerate() {
            if i > 0 {
                DELIMITER.fmt(f)?;
            }

            field_fmt(field, f)?;
        }

        Ok(())
    }
}

/// An error returned when raw attributes fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is invalid.
    InvalidField(field::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField(e) => {
                write!(f, "invalid field")?;

                if let Some(key) = e.key() {
                    write!(f, ": {key}")?;
                }

                Ok(())
            }
        }
    }
}

impl Extend<(Tag, Value)> for Attributes {
    fn extend<T: IntoIterator<Item = (Tag, Value)>>(&mut self, iter: T) {
        self.0.extend(iter);
    }
}

impl FromIterator<(Tag, Value)> for Attributes {
    fn from_iter<T: IntoIterator<Item = (Tag, Value)>>(iter: T) -> Self {
        let mut attributes = Self::default();
        attributes.extend(iter);
        attributes
    }
}

impl FromStr for Attributes {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use self::field::parse_field;

        if s.is_empty() {
            return Ok(Self::default());
        }

        let mut map: IndexMap<Tag, Value> = IndexMap::new();

        for raw_field in s.split(DELIMITER) {
            let (key, value) = parse_field(raw_field).map_err(ParseError::InvalidField)?;

            map.entry(key)
                .and_modify(|v| v.extend(value.iter().cloned()))
                .or_insert(value);
        }

        Ok(Self(map))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let attributes = Attributes::default();
        assert!(attributes.to_string().is_empty());

        let attributes: Attributes = [(Tag::from("gene_id"), Value::from("ndls0"))]
            .into_iter()
            .collect();
        assert_eq!(attributes.to_string(), "gene_id=ndls0");

        let attributes: Attributes = [
            (Tag::from("gene_id"), Value::from("ndls0")),
            (Tag::from("gene_name"), Value::from("gene0")),
        ]
        .into_iter()
        .collect();
        assert_eq!(attributes.to_string(), "gene_id=ndls0;gene_name=gene0")
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let actual: Attributes = "".parse()?;
        let expected = Attributes::default();
        assert_eq!(actual, expected);

        let s = "gene_id=ndls0";
        let actual: Attributes = s.parse()?;
        let expected = [(Tag::from("gene_id"), Value::from("ndls0"))]
            .into_iter()
            .collect();
        assert_eq!(actual, expected);

        let s = "gene_id=ndls0;gene_name=gene0";
        let actual: Attributes = s.parse()?;
        let expected = [
            (Tag::from("gene_id"), Value::from("ndls0")),
            (Tag::from("gene_name"), Value::from("gene0")),
        ]
        .into_iter()
        .collect();
        assert_eq!(actual, expected);

        let s = "gene_id=ndls0;gene_id=ndls1";
        let actual: Attributes = s.parse()?;
        let expected = [(
            Tag::from("gene_id"),
            Value::from(vec![String::from("ndls0"), String::from("ndls1")]),
        )]
        .into_iter()
        .collect();
        assert_eq!(actual, expected);

        Ok(())
    }
}
