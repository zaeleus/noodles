pub mod key;
mod value;

pub use self::{key::Key, value::Value};

use std::{error, fmt};

#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    key: Key,
    value: Value,
}

#[derive(Debug)]
pub struct ParseError(value::ParseError);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid genotype field: {}", self.0)
    }
}

impl Field {
    pub fn from_str_key(s: &str, key: &Key) -> Result<Self, ParseError> {
        Value::from_str_key(s, key)
            .map(|v| Self::new(key.clone(), v))
            .map_err(ParseError)
    }

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
        write!(f, "{}", self.value)
    }
}

#[cfg(test)]
mod tests {
    use crate::header::{format::Type, Number};

    use super::*;

    #[test]
    fn test_from_str_key() -> Result<(), ParseError> {
        let key = Key::ConditionalGenotypeQuality;
        let actual = Field::from_str_key("13", &key)?;
        assert_eq!(actual.key(), &key);
        assert_eq!(actual.value(), &Value::Integer(13));

        let key = Key::Other(String::from("CNQ"), Number::Count(1), Type::Float);
        let actual = Field::from_str_key("8.333", &key)?;
        assert_eq!(actual.key(), &key);
        assert_eq!(actual.value(), &Value::Float(8.333));

        let key = Key::Genotype;
        let actual = Field::from_str_key("0|0", &key)?;
        assert_eq!(actual.key(), &key);
        assert_eq!(actual.value(), &Value::String(String::from("0|0")));

        Ok(())
    }

    #[test]
    fn test_fmt() {
        let field = Field::new(Key::ConditionalGenotypeQuality, Value::Integer(13));
        assert_eq!(field.to_string(), "13");

        let key = Key::Other(String::from("CNQ"), Number::Count(1), Type::Float);
        let field = Field::new(key, Value::Float(8.333));
        assert_eq!(field.to_string(), "8.333");

        let field = Field::new(Key::Genotype, Value::String(String::from("0|0")));
        assert_eq!(field.to_string(), "0|0");
    }
}
