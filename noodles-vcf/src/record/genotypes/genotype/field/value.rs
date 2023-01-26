//! VCF record genotype field value.

pub mod genotype;

pub use self::genotype::Genotype;

use std::{error, fmt, num, str};

use crate::{
    header::{
        format::Type,
        record::value::{map::Format, Map},
        Number,
    },
    record::value::{self, percent_decode},
};

const DELIMITER: char = ',';
const MISSING_VALUE: &str = ".";

/// A VCF record genotype field value.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    /// A 32-bit integer.
    Integer(i32),
    /// A single-precision floating-point.
    Float(f32),
    /// A character.
    Character(char),
    /// A string.
    String(String),
    /// An array of 32-bit integers.
    IntegerArray(Vec<Option<i32>>),
    /// An array of single-precision floating-points.
    FloatArray(Vec<Option<f32>>),
    /// An array of characters.
    CharacterArray(Vec<Option<char>>),
    /// An array of strings.
    StringArray(Vec<Option<String>>),
}

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Integer(n) => write!(f, "{n}"),
            Self::Float(n) => write!(f, "{n}"),
            Self::Character(c) => write!(f, "{c}"),
            Self::String(s) => write!(f, "{s}"),
            Self::IntegerArray(values) => {
                for (i, value) in values.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{DELIMITER}")?;
                    }

                    if let Some(v) = value {
                        write!(f, "{v}")?;
                    } else {
                        f.write_str(MISSING_VALUE)?;
                    }
                }

                Ok(())
            }
            Self::FloatArray(values) => {
                for (i, value) in values.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{DELIMITER}")?;
                    }

                    if let Some(v) = value {
                        write!(f, "{v}")?;
                    } else {
                        f.write_str(MISSING_VALUE)?;
                    }
                }

                Ok(())
            }
            Self::CharacterArray(values) => {
                for (i, value) in values.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{DELIMITER}")?;
                    }

                    if let Some(v) = value {
                        write!(f, "{v}")?;
                    } else {
                        f.write_str(MISSING_VALUE)?;
                    }
                }

                Ok(())
            }
            Self::StringArray(values) => {
                for (i, value) in values.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{DELIMITER}")?;
                    }

                    if let Some(v) = value {
                        write!(f, "{v}")?;
                    } else {
                        f.write_str(MISSING_VALUE)?;
                    }
                }

                Ok(())
            }
        }
    }
}

/// An error returned when a raw VCF record genotype field value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The field cardinality is invalid for the type.
    InvalidNumberForType(Number, Type),
    /// The integer is invalid.
    InvalidInteger(num::ParseIntError),
    /// The floating-point is invalid.
    InvalidFloat(num::ParseFloatError),
    /// The character is invalid.
    InvalidCharacter,
    /// The string is invalid.
    InvalidString(str::Utf8Error),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidInteger(e) => Some(e),
            Self::InvalidFloat(e) => Some(e),
            Self::InvalidString(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidNumberForType(number, ty) => {
                write!(f, "invalid number {number:?} for type {ty:?}")
            }
            Self::InvalidInteger(e) => write!(f, "invalid integer: {e}"),
            Self::InvalidFloat(e) => write!(f, "invalid float: {e}"),
            Self::InvalidCharacter => f.write_str("invalid character"),
            Self::InvalidString(e) => write!(f, "invalid string: {e}"),
        }
    }
}

impl Value {
    /// Parses a raw genotype field value for the given key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::{format::Key, record::value::{map::Format, Map}},
    ///     record::genotypes::genotype::field::Value,
    /// };
    ///
    /// let format = Map::<Format>::from(&Key::ConditionalGenotypeQuality);
    /// assert_eq!(Value::from_str_format("13", &format), Ok(Value::Integer(13)));
    /// ```
    pub fn from_str_format(s: &str, format: &Map<Format>) -> Result<Self, ParseError> {
        match format.ty() {
            Type::Integer => match format.number() {
                Number::Count(0) => Err(ParseError::InvalidNumberForType(
                    format.number(),
                    format.ty(),
                )),
                Number::Count(1) => parse_i32(s),
                _ => parse_i32_array(s),
            },
            Type::Float => match format.number() {
                Number::Count(0) => Err(ParseError::InvalidNumberForType(
                    format.number(),
                    format.ty(),
                )),
                Number::Count(1) => parse_f32(s),
                _ => parse_f32_array(s),
            },
            Type::Character => match format.number() {
                Number::Count(0) => Err(ParseError::InvalidNumberForType(
                    format.number(),
                    format.ty(),
                )),
                Number::Count(1) => parse_char(s),
                _ => parse_char_array(s),
            },
            Type::String => match format.number() {
                Number::Count(0) => Err(ParseError::InvalidNumberForType(
                    format.number(),
                    format.ty(),
                )),
                Number::Count(1) => parse_string(s),
                _ => parse_string_array(s),
            },
        }
    }
}

fn parse_i32(s: &str) -> Result<Value, ParseError> {
    s.parse()
        .map(Value::Integer)
        .map_err(ParseError::InvalidInteger)
}

fn parse_i32_array(s: &str) -> Result<Value, ParseError> {
    s.split(DELIMITER)
        .map(|t| {
            if t == MISSING_VALUE {
                Ok(None)
            } else {
                t.parse().map(Some).map_err(ParseError::InvalidInteger)
            }
        })
        .collect::<Result<_, _>>()
        .map(Value::IntegerArray)
}

fn parse_f32(s: &str) -> Result<Value, ParseError> {
    value::parse_f32(s)
        .map(Value::Float)
        .map_err(ParseError::InvalidFloat)
}

fn parse_f32_array(s: &str) -> Result<Value, ParseError> {
    s.split(DELIMITER)
        .map(|t| {
            if t == MISSING_VALUE {
                Ok(None)
            } else {
                value::parse_f32(t)
                    .map(Some)
                    .map_err(ParseError::InvalidFloat)
            }
        })
        .collect::<Result<_, _>>()
        .map(Value::FloatArray)
}

fn parse_raw_char(s: &str) -> Result<char, ParseError> {
    let mut chars = s.chars();

    if let Some(c) = chars.next() {
        if chars.next().is_none() {
            return Ok(c);
        }
    }

    Err(ParseError::InvalidCharacter)
}

fn parse_char(s: &str) -> Result<Value, ParseError> {
    parse_raw_char(s).map(Value::Character)
}

fn parse_char_array(s: &str) -> Result<Value, ParseError> {
    s.split(DELIMITER)
        .map(|t| {
            if t == MISSING_VALUE {
                Ok(None)
            } else {
                parse_raw_char(t).map(Some)
            }
        })
        .collect::<Result<_, _>>()
        .map(Value::CharacterArray)
}

fn parse_string(s: &str) -> Result<Value, ParseError> {
    percent_decode(s)
        .map(|t| Value::String(t.into()))
        .map_err(ParseError::InvalidString)
}

fn parse_string_array(s: &str) -> Result<Value, ParseError> {
    s.split(DELIMITER)
        .map(|t| {
            if t == MISSING_VALUE {
                Ok(None)
            } else {
                percent_decode(t)
                    .map(|u| Some(u.into()))
                    .map_err(ParseError::InvalidString)
            }
        })
        .collect::<Result<_, _>>()
        .map(Value::StringArray)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let value = Value::Integer(2);
        assert_eq!(value.to_string(), "2");

        let value = Value::Float(0.333);
        assert_eq!(value.to_string(), "0.333");

        let value = Value::Character('n');
        assert_eq!(value.to_string(), "n");

        let value = Value::String(String::from("noodles"));
        assert_eq!(value.to_string(), "noodles");

        let value = Value::IntegerArray(vec![Some(2)]);
        assert_eq!(value.to_string(), "2");

        let value = Value::IntegerArray(vec![Some(2), Some(5)]);
        assert_eq!(value.to_string(), "2,5");

        let value = Value::IntegerArray(vec![Some(2), None]);
        assert_eq!(value.to_string(), "2,.");

        let value = Value::FloatArray(vec![Some(0.333)]);
        assert_eq!(value.to_string(), "0.333");

        let value = Value::FloatArray(vec![Some(0.333), Some(0.667)]);
        assert_eq!(value.to_string(), "0.333,0.667");

        let value = Value::FloatArray(vec![Some(0.333), None]);
        assert_eq!(value.to_string(), "0.333,.");

        let value = Value::CharacterArray(vec![Some('n')]);
        assert_eq!(value.to_string(), "n");

        let value = Value::CharacterArray(vec![Some('n'), Some('d'), Some('l'), Some('s')]);
        assert_eq!(value.to_string(), "n,d,l,s");

        let value = Value::CharacterArray(vec![Some('n'), Some('d'), Some('l'), None]);
        assert_eq!(value.to_string(), "n,d,l,.");

        let value = Value::StringArray(vec![Some(String::from("noodles"))]);
        assert_eq!(value.to_string(), "noodles");

        let value = Value::StringArray(vec![
            Some(String::from("noodles")),
            Some(String::from("vcf")),
        ]);
        assert_eq!(value.to_string(), "noodles,vcf");

        let value = Value::StringArray(vec![Some(String::from("noodles")), None]);
        assert_eq!(value.to_string(), "noodles,.");
    }

    #[test]
    fn test_from_str_format_with_integer() -> Result<(), crate::header::format::key::ParseError> {
        let format = Map::<Format>::new(Number::Count(0), Type::Integer, String::new());
        assert_eq!(
            Value::from_str_format("8", &format),
            Err(ParseError::InvalidNumberForType(
                Number::Count(0),
                Type::Integer
            ))
        );

        let format = Map::<Format>::new(Number::Count(1), Type::Integer, String::new());
        assert_eq!(Value::from_str_format("8", &format), Ok(Value::Integer(8)));

        let format = Map::<Format>::new(Number::Count(2), Type::Integer, String::new());
        assert_eq!(
            Value::from_str_format("8,13", &format),
            Ok(Value::IntegerArray(vec![Some(8), Some(13)]))
        );

        let format = Map::<Format>::new(Number::Count(2), Type::Integer, String::new());
        assert_eq!(
            Value::from_str_format("8,.", &format),
            Ok(Value::IntegerArray(vec![Some(8), None]))
        );

        Ok(())
    }

    #[test]
    fn test_from_str_format_with_float() -> Result<(), crate::header::format::key::ParseError> {
        let format = Map::<Format>::new(Number::Count(0), Type::Float, String::new());
        assert_eq!(
            Value::from_str_format("0.333", &format),
            Err(ParseError::InvalidNumberForType(
                Number::Count(0),
                Type::Float
            ))
        );

        let format = Map::<Format>::new(Number::Count(1), Type::Float, String::new());
        assert_eq!(
            Value::from_str_format("0.333", &format),
            Ok(Value::Float(0.333))
        );

        let format = Map::<Format>::new(Number::Count(2), Type::Float, String::new());
        assert_eq!(
            Value::from_str_format("0.333,0.667", &format),
            Ok(Value::FloatArray(vec![Some(0.333), Some(0.667)]))
        );

        let format = Map::<Format>::new(Number::Count(2), Type::Float, String::new());
        assert_eq!(
            Value::from_str_format("0.333,.", &format),
            Ok(Value::FloatArray(vec![Some(0.333), None]))
        );

        Ok(())
    }

    #[test]
    fn test_from_str_format_with_character() -> Result<(), crate::header::format::key::ParseError> {
        let format = Map::<Format>::new(Number::Count(0), Type::Character, String::new());
        assert_eq!(
            Value::from_str_format("n", &format),
            Err(ParseError::InvalidNumberForType(
                Number::Count(0),
                Type::Character
            ))
        );

        let format = Map::<Format>::new(Number::Count(1), Type::Character, String::new());
        assert_eq!(
            Value::from_str_format("n", &format),
            Ok(Value::Character('n'))
        );

        let format = Map::<Format>::new(Number::Count(2), Type::Character, String::new());
        assert_eq!(
            Value::from_str_format("n,d,l,s", &format),
            Ok(Value::CharacterArray(vec![
                Some('n'),
                Some('d'),
                Some('l'),
                Some('s')
            ]))
        );

        let format = Map::<Format>::new(Number::Count(2), Type::Character, String::new());
        assert_eq!(
            Value::from_str_format("n,d,l,.", &format),
            Ok(Value::CharacterArray(vec![
                Some('n'),
                Some('d'),
                Some('l'),
                None
            ]))
        );

        Ok(())
    }

    #[test]
    fn test_from_str_format_with_string() -> Result<(), crate::header::format::key::ParseError> {
        let format = Map::<Format>::new(Number::Count(0), Type::String, String::new());
        assert_eq!(
            Value::from_str_format("noodles", &format),
            Err(ParseError::InvalidNumberForType(
                Number::Count(0),
                Type::String
            ))
        );

        let format = Map::<Format>::new(Number::Count(1), Type::String, String::new());
        assert_eq!(
            Value::from_str_format("noodles", &format),
            Ok(Value::String(String::from("noodles")))
        );
        assert_eq!(
            Value::from_str_format("8%25", &format),
            Ok(Value::String(String::from("8%")))
        );

        let format = Map::<Format>::new(Number::Count(2), Type::String, String::new());
        assert_eq!(
            Value::from_str_format("noodles,vcf", &format),
            Ok(Value::StringArray(vec![
                Some(String::from("noodles")),
                Some(String::from("vcf"))
            ]))
        );
        assert_eq!(
            Value::from_str_format("8%25,13%25", &format),
            Ok(Value::StringArray(vec![
                Some(String::from("8%")),
                Some(String::from("13%")),
            ]))
        );
        assert_eq!(
            Value::from_str_format("noodles,.", &format),
            Ok(Value::StringArray(vec![
                Some(String::from("noodles")),
                None,
            ]))
        );

        Ok(())
    }
}
