//! VCF header record and components.

pub mod key;
pub(crate) mod parser;
mod value;

use indexmap::IndexMap;

pub use self::{key::Key, value::Value};

use std::{error, fmt, str::FromStr};

pub(crate) const PREFIX: &str = "##";

/// A generic VCF header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Record {
    key: Key,
    value: Value,
}

impl Record {
    /// Creates a generic VCF header record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{record::{key, Value}, Record};
    /// let record = Record::new(key::FILE_FORMAT, Value::from("VCFv4.3"));
    /// ```
    pub fn new(key: Key, value: Value) -> Self {
        Self { key, value }
    }

    /// Returns the key of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{record::{key, Value}, Record};
    /// let record = Record::new(key::FILE_FORMAT, Value::from("VCFv4.3"));
    /// assert_eq!(record.key(), &key::FILE_FORMAT);
    /// ```
    pub fn key(&self) -> &Key {
        &self.key
    }

    /// Returns the value of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{record::{key, Value}, Record};
    /// let record = Record::new(key::FILE_FORMAT, Value::from("VCFv4.3"));
    /// assert_eq!(record.value(), &Value::from("VCFv4.3"));
    /// ```
    pub fn value(&self) -> &Value {
        &self.value
    }
}

/// An error returned when a raw VCF header record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (_, (raw_key, raw_value)) = parser::parse(s).map_err(|_| ParseError::Invalid)?;

        let key = Key::from(raw_key);

        let value = match raw_value {
            parser::Value::String(s) => Value::String(s),
            parser::Value::Struct(raw_fields) => {
                let mut fields = IndexMap::new();

                for (k, v) in raw_fields {
                    if fields.insert(k, v).is_some() {
                        return Err(ParseError::Invalid);
                    }
                }

                let id = match fields.shift_remove("ID") {
                    Some(value) => value,
                    None => return Err(ParseError::Invalid),
                };

                Value::Struct(id, fields)
            }
        };

        Ok(Self::new(key, value))
    }
}

impl From<Record> for (Key, Value) {
    fn from(record: Record) -> Self {
        (record.key, record.value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        let line = "##fileformat=VCFv4.3";

        assert_eq!(
            line.parse(),
            Ok(Record::new(key::FILE_FORMAT, Value::from("VCFv4.3")))
        );

        let line =
            r#"##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">"#;

        assert_eq!(
            line.parse(),
            Ok(Record::new(
                key::INFO,
                Value::Struct(
                    String::from("NS"),
                    [
                        (String::from("Number"), String::from("1")),
                        (String::from("Type"), String::from("Integer")),
                        (
                            String::from("Description"),
                            String::from("Number of samples with data"),
                        ),
                    ]
                    .into_iter()
                    .collect()
                )
            ))
        );

        assert!("".parse::<Record>().is_err());
    }
}
