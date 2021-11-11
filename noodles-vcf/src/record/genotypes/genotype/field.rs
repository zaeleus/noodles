//! VCF record genotype field.

pub mod key;
pub mod value;

pub use self::{key::Key, value::Value};

use std::{error, fmt};

const MISSING_VALUE: &str = ".";

/// A VCF record genotype field.
#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    key: Key,
    value: Option<Value>,
}

/// An error returned when a raw VCF record info field value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The value is invalid.
    InvalidValue(value::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidValue(e) => write!(f, "invalid value: {}", e),
        }
    }
}

impl Field {
    /// Parses a raw genotype field for the given key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::genotypes::genotype::{field::{Key, Value}, Field};
    ///
    /// assert_eq!(
    ///     Field::from_str_key("13", &Key::ConditionalGenotypeQuality),
    ///     Ok(Field::new(Key::ConditionalGenotypeQuality, Some(Value::Integer(13))))
    /// );
    /// ```
    pub fn from_str_key(s: &str, key: &Key) -> Result<Self, ParseError> {
        if s == MISSING_VALUE {
            Ok(Self::new(key.clone(), None))
        } else {
            Value::from_str_key(s, key)
                .map(|v| Self::new(key.clone(), Some(v)))
                .map_err(ParseError::InvalidValue)
        }
    }

    /// Creates a VCF record genotype field.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::genotypes::genotype::{field::{Key, Value}, Field};
    ///
    /// let field = Field::new(
    ///     Key::ConditionalGenotypeQuality,
    ///     Some(Value::Integer(13)),
    /// );
    /// ```
    pub fn new(key: Key, value: Option<Value>) -> Self {
        Self { key, value }
    }

    /// Returns the genotype field key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::genotypes::genotype::{field::{Key, Value}, Field};
    ///
    /// let field = Field::new(
    ///     Key::ConditionalGenotypeQuality,
    ///     Some(Value::Integer(13)),
    /// );
    ///
    /// assert_eq!(field.key(), &Key::ConditionalGenotypeQuality);
    /// ```
    pub fn key(&self) -> &Key {
        &self.key
    }

    /// Returns the genotype field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::genotypes::genotype::{field::{Key, Value}, Field};
    ///
    /// let field = Field::new(
    ///     Key::ConditionalGenotypeQuality,
    ///     Some(Value::Integer(13)),
    /// );
    ///
    /// assert_eq!(field.value(), Some(&Value::Integer(13)));
    /// ```
    pub fn value(&self) -> Option<&Value> {
        self.value.as_ref()
    }

    /// Returns a mutable reference to the value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::genotypes::genotype::{field::{Key, Value}, Field};
    /// let mut field = Field::new(Key::ConditionalGenotypeQuality, Some(Value::Integer(13)));
    /// *field.value_mut() = None;
    /// assert!(field.value().is_none());
    /// ```
    pub fn value_mut(&mut self) -> &mut Option<Value> {
        &mut self.value
    }
}

impl fmt::Display for Field {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(value) = self.value() {
            write!(f, "{}", value)
        } else {
            f.write_str(MISSING_VALUE)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::header::{format::Type, Number};

    use super::*;

    #[test]
    fn test_from_str_key() -> Result<(), ParseError> {
        let key = Key::MappingQuality;
        let actual = Field::from_str_key(".", &key)?;
        assert_eq!(actual.key(), &key);
        assert_eq!(actual.value(), None);

        let key = Key::ConditionalGenotypeQuality;
        let actual = Field::from_str_key("13", &key)?;
        assert_eq!(actual.key(), &key);
        assert_eq!(actual.value(), Some(&Value::Integer(13)));

        let key = Key::Other(
            String::from("CNQ"),
            Number::Count(1),
            Type::Float,
            String::default(),
        );
        let actual = Field::from_str_key("8.333", &key)?;
        assert_eq!(actual.key(), &key);
        assert_eq!(actual.value(), Some(&Value::Float(8.333)));

        let key = Key::Genotype;
        let actual = Field::from_str_key("0|0", &key)?;
        assert_eq!(actual.key(), &key);
        assert_eq!(actual.value(), Some(&Value::String(String::from("0|0"))));

        Ok(())
    }

    #[test]
    fn test_fmt() {
        let field = Field::new(Key::MappingQuality, None);
        assert_eq!(field.to_string(), ".");

        let field = Field::new(Key::ConditionalGenotypeQuality, Some(Value::Integer(13)));
        assert_eq!(field.to_string(), "13");

        let key = Key::Other(
            String::from("CNQ"),
            Number::Count(1),
            Type::Float,
            String::default(),
        );
        let field = Field::new(key, Some(Value::Float(8.333)));
        assert_eq!(field.to_string(), "8.333");

        let field = Field::new(Key::Genotype, Some(Value::String(String::from("0|0"))));
        assert_eq!(field.to_string(), "0|0");
    }
}
