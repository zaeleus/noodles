//! SAM record data and fields.

pub mod field;

pub use self::field::Field;

use std::{
    convert::TryFrom,
    error,
    fmt::{self, Write},
    ops::{Deref, DerefMut},
    str::FromStr,
};

use indexmap::IndexMap;

const DELIMITER: char = '\t';

/// SAM record data.
///
/// This is also called optional fields.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Data(IndexMap<field::Tag, Field>);

impl Deref for Data {
    type Target = IndexMap<field::Tag, Field>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Data {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl fmt::Display for Data {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, field) in self.values().enumerate() {
            if i > 0 {
                f.write_char('\t')?;
            }
            write!(f, "{}", field)?;
        }

        Ok(())
    }
}

/// An error returned when raw SAM record data fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input data contains an invalid field.
    InvalidField(field::ParseError),
    /// A tag is duplicated.
    ///
    /// § 1.5 The alignment section: optional fields (2021-01-07): "Each `TAG` can only appear once
    /// in one alignment line."
    DuplicateTag(field::Tag),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField(e) => write!(f, "invalid field: {}", e),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {}", tag),
        }
    }
}

impl FromStr for Data {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Ok(Self::default());
        }

        let mut data = Self::default();

        for s in s.split(DELIMITER) {
            let field: Field = s.parse().map_err(ParseError::InvalidField)?;
            let tag = field.tag();

            if data.insert(tag, field).is_some() {
                return Err(ParseError::DuplicateTag(tag));
            }
        }

        Ok(data)
    }
}

impl TryFrom<Vec<Field>> for Data {
    type Error = ParseError;

    fn try_from(fields: Vec<Field>) -> Result<Self, Self::Error> {
        let mut map = IndexMap::new();

        for field in fields {
            let tag = field.tag();

            if map.insert(tag, field).is_some() {
                return Err(ParseError::DuplicateTag(tag));
            }
        }

        Ok(Self(map))
    }
}

#[cfg(test)]
mod tests {
    use super::field::{Tag, Value};

    use super::*;

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let data = Data::try_from(vec![
            Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
            Field::new(Tag::AlignmentHitCount, Value::Int(1)),
        ])?;

        let expected = "RG:Z:rg0\tNH:i:1";

        assert_eq!(data.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("".parse(), Ok(Data::default()));

        assert_eq!(
            "RG:Z:rg0\tNH:i:1".parse(),
            Ok(Data::try_from(vec![
                Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
                Field::new(Tag::AlignmentHitCount, Value::Int(1)),
            ])?)
        );

        assert_eq!(
            "NH:i:1\tNH:i:1".parse::<Data>(),
            Err(ParseError::DuplicateTag(Tag::AlignmentHitCount))
        );

        Ok(())
    }
}
