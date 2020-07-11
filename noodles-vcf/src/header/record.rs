pub mod kind;
mod parser;
mod value;

pub use self::{kind::Kind, value::Value};

use std::{error, fmt, str::FromStr};

/// A generic VCF header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Record {
    key: Kind,
    value: Value,
}

impl Record {
    /// Creates a generic VCF header record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{record::{Kind, Value}, Record};
    /// let record = Record::new(Kind::FileFormat, Value::String(String::from("VCFv4.3")));
    /// ```
    pub fn new(key: Kind, value: Value) -> Self {
        Self { key, value }
    }

    /// Returns the key of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{record::{Kind, Value}, Record};
    /// let record = Record::new(Kind::FileFormat, Value::String(String::from("VCFv4.3")));
    /// assert_eq!(record.key(), &Kind::FileFormat);
    /// ```
    pub fn key(&self) -> &Kind {
        &self.key
    }

    /// Returns the value of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{record::{Kind, Value}, Record};
    /// let record = Record::new(Kind::FileFormat, Value::String(String::from("VCFv4.3")));
    /// assert_eq!(record.value(), &Value::String(String::from("VCFv4.3")));
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
    /// The record key is invalid.
    InvalidKind(kind::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid record: ")?;

        match self {
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidKind(e) => write!(f, "invalid kind: {}", e),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (_, (key, value)) = parser::parse(s).map_err(|_| ParseError::Invalid)?;
        let kind = key.parse().map_err(ParseError::InvalidKind)?;
        Ok(Self::new(kind, value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        let line =
            r#"##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">"#;

        assert_eq!(
            line.parse(),
            Ok(Record::new(
                Kind::Info,
                Value::Struct(vec![
                    (String::from("ID"), String::from("NS")),
                    (String::from("Number"), String::from("1")),
                    (String::from("Type"), String::from("Integer")),
                    (
                        String::from("Description"),
                        String::from("Number of samples with data"),
                    ),
                ])
            ))
        );

        assert!("".parse::<Record>().is_err());
    }
}
