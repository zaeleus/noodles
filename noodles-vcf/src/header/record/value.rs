//! VCF header record value.

use std::{
    error,
    fmt::{self, Write},
};

use indexmap::IndexMap;

/// A VCF header record value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Value {
    /// A string.
    String(String),
    /// A structure.
    Struct(String, IndexMap<String, String>),
}

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::String(value) => f.write_str(value),
            Self::Struct(id, fields) => {
                f.write_char('<')?;

                write!(f, "ID={}", id)?;

                for (key, value) in fields {
                    write!(f, ",{}={}", key, value)?;
                }

                f.write_char('>')?;

                Ok(())
            }
        }
    }
}

impl From<&str> for Value {
    fn from(s: &str) -> Self {
        Self::String(s.into())
    }
}

impl From<String> for Value {
    fn from(s: String) -> Self {
        Self::String(s)
    }
}

/// An error returned when raw fields fail to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromFieldsError {
    /// The ID field is missing.
    MissingId,
    /// A tag is duplicated.
    DuplicateTag,
}

impl error::Error for TryFromFieldsError {}

impl fmt::Display for TryFromFieldsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingId => "missing ID".fmt(f),
            Self::DuplicateTag => "duplicate tag".fmt(f),
        }
    }
}

impl TryFrom<Vec<(String, String)>> for Value {
    type Error = TryFromFieldsError;

    fn try_from(fields: Vec<(String, String)>) -> Result<Self, Self::Error> {
        let mut map = IndexMap::new();

        for (key, value) in fields {
            if map.insert(key, value).is_some() {
                return Err(TryFromFieldsError::DuplicateTag);
            }
        }

        let id = map
            .shift_remove("ID")
            .ok_or(TryFromFieldsError::MissingId)?;

        Ok(Self::Struct(id, map))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Value::from("VCFv4.3").to_string(), "VCFv4.3");

        assert_eq!(
            Value::Struct(String::from("sq0"), IndexMap::new()).to_string(),
            "<ID=sq0>"
        );

        assert_eq!(
            Value::Struct(
                String::from("sq0"),
                [(String::from("length"), String::from("13"))]
                    .into_iter()
                    .collect(),
            )
            .to_string(),
            "<ID=sq0,length=13>"
        );
    }

    #[test]
    fn test_from_str_for_value() {
        assert_eq!(
            Value::from("VCFv4.3"),
            Value::String(String::from("VCFv4.3"))
        );
    }

    #[test]
    fn test_from_string_for_value() {
        assert_eq!(
            Value::from(String::from("VCFv4.3")),
            Value::String(String::from("VCFv4.3"))
        );
    }

    #[test]
    fn test_from_vec_string_string_for_value() {
        assert_eq!(
            Value::try_from(vec![(String::from("ID"), String::from("sq0"))]),
            Ok(Value::Struct(String::from("sq0"), Default::default(),))
        );

        assert_eq!(
            Value::try_from(vec![
                (String::from("ID"), String::from("sq0")),
                (String::from("length"), String::from("13"))
            ]),
            Ok(Value::Struct(
                String::from("sq0"),
                [(String::from("length"), String::from("13"))]
                    .into_iter()
                    .collect(),
            ))
        );

        assert_eq!(
            Value::try_from(vec![(String::from("length"), String::from("13"))]),
            Err(TryFromFieldsError::MissingId),
        );

        assert_eq!(
            Value::try_from(vec![
                (String::from("ID"), String::from("sq0")),
                (String::from("ID"), String::from("sq0")),
            ]),
            Err(TryFromFieldsError::DuplicateTag),
        );
    }
}
