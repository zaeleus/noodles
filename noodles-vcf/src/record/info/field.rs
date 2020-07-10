pub mod key;
pub mod value;

pub use self::{key::Key, value::Value};

use std::{error, fmt, str::FromStr};

use crate::header::info::Type;

const SEPARATOR: char = '=';
const MAX_COMPONENTS: usize = 2;

#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    key: Key,
    value: Value,
}

impl Field {
    pub fn new(key: Key, value: Value) -> Self {
        Self { key, value }
    }

    pub fn key(&self) -> &Key {
        &self.key
    }

    pub fn value(&self) -> &Value {
        &self.value
    }
}

impl fmt::Display for Field {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.value() {
            Value::Flag => write!(f, "{}", self.key),
            value => write!(f, "{}{}{}", self.key, SEPARATOR, value),
        }
    }
}

/// An error returned when a raw VCF record info field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The key is missing.
    MissingKey,
    /// The key is invalid.
    InvalidKey(key::ParseError),
    /// The value is missing.
    MissingValue,
    /// The value is invalid.
    InvalidValue(value::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid info: ")?;

        match self {
            ParseError::MissingKey => f.write_str("missing key"),
            ParseError::InvalidKey(e) => write!(f, "invalid key: {}", e),
            ParseError::MissingValue => f.write_str("missing value"),
            ParseError::InvalidValue(e) => write!(f, "invalid value: {}", e),
        }
    }
}

impl FromStr for Field {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut components = s.splitn(MAX_COMPONENTS, SEPARATOR);

        let key: Key = components
            .next()
            .ok_or_else(|| ParseError::MissingKey)
            .and_then(|s| s.parse().map_err(ParseError::InvalidKey))?;

        let value = if let Type::Flag = key.ty() {
            let s = components.next().unwrap_or_default();
            Value::from_str_key(s, &key).map_err(ParseError::InvalidValue)?
        } else if let Key::Other(..) = key {
            if let Some(s) = components.next() {
                Value::from_str_key(s, &key).map_err(ParseError::InvalidValue)?
            } else {
                Value::Flag
            }
        } else {
            components
                .next()
                .ok_or_else(|| ParseError::MissingValue)
                .and_then(|s| Value::from_str_key(s, &key).map_err(ParseError::InvalidValue))?
        };

        Ok(Self::new(key, value))
    }
}

#[cfg(test)]
mod tests {
    use crate::header::Number;

    use super::*;

    #[test]
    fn test_fmt() {
        let field = Field::new(Key::SamplesWithDataCount, Value::Integer(2));
        assert_eq!(field.to_string(), "NS=2");

        let field = Field::new(Key::BaseQuality, Value::Float(1.333));
        assert_eq!(field.to_string(), "BQ=1.333");

        let field = Field::new(Key::IsSomaticMutation, Value::Flag);
        assert_eq!(field.to_string(), "SOMATIC");

        let field = Field::new(
            Key::Other(String::from("SVTYPE"), Number::Count(1), Type::String),
            Value::String(String::from("DEL")),
        );
        assert_eq!(field.to_string(), "SVTYPE=DEL");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let actual: Field = "NS=2".parse()?;
        assert_eq!(actual.key(), &Key::SamplesWithDataCount);
        assert_eq!(actual.value(), &Value::Integer(2));

        let actual: Field = "BQ=1.333".parse()?;
        assert_eq!(actual.key(), &Key::BaseQuality);
        assert_eq!(actual.value(), &Value::Float(1.333));

        let actual: Field = "SOMATIC".parse()?;
        assert_eq!(actual.key(), &Key::IsSomaticMutation);
        assert_eq!(actual.value(), &Value::Flag);

        let actual: Field = "EVENT=INV0".parse()?;
        assert_eq!(actual.key(), &Key::BreakendEventId);
        assert_eq!(actual.value(), &Value::String(String::from("INV0")));

        let actual: Field = "NDLS=VCF".parse()?;
        assert_eq!(
            actual.key(),
            &Key::Other(String::from("NDLS"), Number::Count(1), Type::String)
        );
        assert_eq!(actual.value(), &Value::String(String::from("VCF")));

        let actual: Field = "FLG".parse()?;
        assert_eq!(
            actual.key(),
            &Key::Other(String::from("FLG"), Number::Count(1), Type::String)
        );
        assert_eq!(actual.value(), &Value::Flag);

        Ok(())
    }
}
