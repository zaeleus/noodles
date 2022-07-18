//! VCF record info field.

pub mod value;

pub use self::value::Value;

use std::{error, fmt, str::FromStr};

use crate::header::{
    self,
    info::{Key, Type},
    record::value::{map::Info, Map},
    Infos,
};

const MISSING_VALUE: &str = ".";
const SEPARATOR: char = '=';

/// A VCF record info field.
#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    key: Key,
    value: Option<Value>,
}

impl Field {
    /// Parses a raw VCF record info field.
    pub fn try_from_str(s: &str, infos: &Infos) -> Result<Self, ParseError> {
        parse(s, infos)
    }

    /// Creates a VCF record info field.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::Key,
    ///     record::info::{field::Value, Field},
    /// };
    ///
    /// let field = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(1)));
    /// ```
    pub fn new(key: Key, value: Option<Value>) -> Self {
        Self { key, value }
    }

    /// Returns the field key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::Key,
    ///     record::info::{field::Value, Field},
    /// };
    ///
    /// let field = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(1)));
    /// assert_eq!(field.key(), &Key::SamplesWithDataCount);
    /// ```
    pub fn key(&self) -> &Key {
        &self.key
    }

    /// Returns the field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::Key,
    ///     record::info::{field::Value, Field},
    /// };
    ///
    /// let field = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(1)));
    /// assert_eq!(field.value(), Some(&Value::Integer(1)));
    /// ```
    pub fn value(&self) -> Option<&Value> {
        self.value.as_ref()
    }

    /// Returns a mutable reference to the value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::Key,
    ///     record::info::{field::Value, Field},
    /// };
    ///
    /// let mut field = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(1)));
    /// *field.value_mut() = Some(Value::Integer(2));
    /// assert_eq!(field.value(), Some(&Value::Integer(2)));
    /// ```
    pub fn value_mut(&mut self) -> &mut Option<Value> {
        &mut self.value
    }
}

impl fmt::Display for Field {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.value() {
            None => write!(f, "{}{}{}", self.key, SEPARATOR, MISSING_VALUE),
            Some(Value::Flag) => write!(f, "{}", self.key),
            Some(value) => write!(f, "{}{}{}", self.key, SEPARATOR, value),
        }
    }
}

/// An error returned when a raw VCF record info field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The key is missing.
    MissingKey,
    /// The key is invalid.
    InvalidKey(header::info::key::ParseError),
    /// The value is missing.
    MissingValue,
    /// The value is invalid.
    InvalidValue(value::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingKey => f.write_str("missing key"),
            Self::InvalidKey(e) => write!(f, "invalid key: {}", e),
            Self::MissingValue => f.write_str("missing value"),
            Self::InvalidValue(e) => write!(f, "invalid value: {}", e),
        }
    }
}

impl FromStr for Field {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::try_from_str(s, &Infos::default())
    }
}

fn parse(s: &str, infos: &Infos) -> Result<Field, ParseError> {
    const MAX_COMPONENTS: usize = 2;

    let mut components = s.splitn(MAX_COMPONENTS, SEPARATOR);

    let key: Key = components
        .next()
        .ok_or(ParseError::MissingKey)
        .and_then(|t| t.parse().map_err(ParseError::InvalidKey))?;

    let value = if let Some(info) = infos.get(&key) {
        parse_value(&mut components, info)?
    } else {
        let info = Map::<Info>::from(key.clone());
        parse_value(&mut components, &info)?
    };

    Ok(Field::new(key, value))
}

fn parse_value<'a, I>(iter: &mut I, info: &Map<Info>) -> Result<Option<Value>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    if let Type::Flag = info.ty() {
        let t = iter.next().unwrap_or_default();

        if t == MISSING_VALUE {
            Ok(None)
        } else {
            Value::from_str_info(t, info)
                .map(Some)
                .map_err(ParseError::InvalidValue)
        }
    } else if let Key::Other(..) = info.id() {
        if let Some(t) = iter.next() {
            if t == MISSING_VALUE {
                Ok(None)
            } else {
                Value::from_str_info(t, info)
                    .map(Some)
                    .map_err(ParseError::InvalidValue)
            }
        } else {
            Ok(Some(Value::Flag))
        }
    } else if let Some(t) = iter.next() {
        if t == MISSING_VALUE {
            Ok(None)
        } else {
            Value::from_str_info(t, info)
                .map(Some)
                .map_err(ParseError::InvalidValue)
        }
    } else {
        Err(ParseError::MissingValue)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() -> Result<(), crate::header::info::key::ParseError> {
        let field = Field::new(Key::AlleleCount, None);
        assert_eq!(field.to_string(), "AC=.");

        let field = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
        assert_eq!(field.to_string(), "NS=2");

        let field = Field::new(Key::BaseQuality, Some(Value::Float(1.333)));
        assert_eq!(field.to_string(), "BQ=1.333");

        let field = Field::new(Key::IsSomaticMutation, Some(Value::Flag));
        assert_eq!(field.to_string(), "SOMATIC");

        let field = Field::new("NOODLES".parse()?, Some(Value::String(String::from("VCF"))));
        assert_eq!(field.to_string(), "NOODLES=VCF");

        Ok(())
    }

    #[test]
    fn test_parse() -> Result<(), crate::header::info::key::ParseError> {
        let header = crate::Header::builder()
            .add_info(Map::<Info>::from(Key::AlleleCount))
            .add_info(Map::<Info>::from(Key::SamplesWithDataCount))
            .add_info(Map::<Info>::from(Key::IsSomaticMutation))
            .add_info(Map::<Info>::from(Key::BreakendEventId))
            .build();

        assert_eq!(
            parse("AC=.", header.infos()),
            Ok(Field::new(Key::AlleleCount, None))
        );

        assert_eq!(
            parse("NS=2", header.infos()),
            Ok(Field::new(
                Key::SamplesWithDataCount,
                Some(Value::Integer(2))
            ))
        );

        assert_eq!(
            parse("BQ=1.333", header.infos()),
            Ok(Field::new(Key::BaseQuality, Some(Value::Float(1.333))))
        );

        assert_eq!(
            parse("SOMATIC", header.infos()),
            Ok(Field::new(Key::IsSomaticMutation, Some(Value::Flag)))
        );

        assert_eq!(
            parse("EVENT=INV0", header.infos()),
            Ok(Field::new(
                Key::BreakendEventId,
                Some(Value::String(String::from("INV0")))
            ))
        );

        let key = "NDLS".parse()?;
        assert_eq!(
            parse("NDLS=VCF", header.infos()),
            Ok(Field::new(key, Some(Value::String(String::from("VCF")))))
        );

        let key = "FLG".parse()?;
        assert_eq!(
            parse("FLG", header.infos()),
            Ok(Field::new(key, Some(Value::Flag)))
        );

        Ok(())
    }
}
