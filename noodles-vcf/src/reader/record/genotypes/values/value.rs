use std::{error, fmt, num, str};

use noodles_core as core;

use crate::{
    header::{record::value::map::format::Type, Number},
    reader::record::MISSING,
    record::{
        genotypes::sample::{value::Array, Value},
        value,
    },
};

const DELIMITER: char = ',';

/// An error returned when a raw VCF record genotypes values value fails to parse.
#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The field cardinality is invalid for the type.
    InvalidNumberForType(Number, Type),
    /// The integer value is invalid.
    InvalidInteger(num::ParseIntError),
    /// The float value is invalid.
    InvalidFloat(num::ParseFloatError),
    /// The character value is invalid.
    InvalidCharacter,
    /// The string value is invalid.
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
            ParseError::InvalidNumberForType(number, ty) => {
                write!(f, "invalid number {number:?} for type {ty:?}")
            }
            ParseError::InvalidInteger(_) => write!(f, "invalid integer"),
            ParseError::InvalidFloat(_) => write!(f, "invalid float"),
            ParseError::InvalidCharacter => write!(f, "invalid character"),
            ParseError::InvalidString(_) => write!(f, "invalid string"),
        }
    }
}

impl From<ParseError> for core::Error {
    fn from(e: ParseError) -> Self {
        Self::new(core::error::Kind::Parse, e)
    }
}

pub(super) fn parse_value(number: Number, ty: Type, s: &str) -> Result<Value, ParseError> {
    match (number, ty) {
        (Number::Count(0), _) => Err(ParseError::InvalidNumberForType(number, ty)),
        (Number::Count(1), Type::Integer) => parse_i32(s),
        (Number::Count(1), Type::Float) => parse_f32(s),
        (Number::Count(1), Type::Character) => parse_char(s),
        (Number::Count(1), Type::String) => parse_string(s),
        (_, Type::Integer) => parse_i32_array(s),
        (_, Type::Float) => parse_f32_array(s),
        (_, Type::Character) => parse_char_array(s),
        (_, Type::String) => parse_string_array(s),
    }
}

fn parse_i32(s: &str) -> Result<Value, ParseError> {
    s.parse()
        .map(Value::Integer)
        .map_err(ParseError::InvalidInteger)
}

fn parse_i32_array(s: &str) -> Result<Value, ParseError> {
    s.split(DELIMITER)
        .map(|t| match t {
            MISSING => Ok(None),
            _ => t.parse().map(Some).map_err(ParseError::InvalidInteger),
        })
        .collect::<Result<_, _>>()
        .map(|values| Value::Array(Array::Integer(values)))
}

fn parse_f32(s: &str) -> Result<Value, ParseError> {
    value::parse_f32(s)
        .map(Value::Float)
        .map_err(ParseError::InvalidFloat)
}

fn parse_f32_array(s: &str) -> Result<Value, ParseError> {
    s.split(DELIMITER)
        .map(|t| match t {
            MISSING => Ok(None),
            _ => t.parse().map(Some).map_err(ParseError::InvalidFloat),
        })
        .collect::<Result<_, _>>()
        .map(|values| Value::Array(Array::Float(values)))
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
        .map(|t| match t {
            MISSING => Ok(None),
            _ => parse_raw_char(t).map(Some),
        })
        .collect::<Result<_, _>>()
        .map(|values| Value::Array(Array::Character(values)))
}

fn parse_raw_string(s: &str) -> Result<String, ParseError> {
    value::percent_decode(s)
        .map(|t| t.into())
        .map_err(ParseError::InvalidString)
}

fn parse_string(s: &str) -> Result<Value, ParseError> {
    parse_raw_string(s).map(Value::String)
}

fn parse_string_array(s: &str) -> Result<Value, ParseError> {
    s.split(DELIMITER)
        .map(|t| match t {
            MISSING => Ok(None),
            _ => parse_raw_string(t).map(Some),
        })
        .collect::<Result<_, _>>()
        .map(|values| Value::Array(Array::String(values)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_value_with_integer() {
        assert_eq!(
            parse_value(Number::Count(0), Type::Integer, "8"),
            Err(ParseError::InvalidNumberForType(
                Number::Count(0),
                Type::Integer
            ))
        );

        assert_eq!(
            parse_value(Number::Count(1), Type::Integer, "8"),
            Ok(Value::Integer(8))
        );

        assert_eq!(
            parse_value(Number::Count(2), Type::Integer, "8,13"),
            Ok(Value::Array(Array::Integer(vec![Some(8), Some(13)])))
        );

        assert_eq!(
            parse_value(Number::Count(2), Type::Integer, "8,."),
            Ok(Value::Array(Array::Integer(vec![Some(8), None])))
        );
    }

    #[test]
    fn test_parse_value_with_float() {
        assert_eq!(
            parse_value(Number::Count(0), Type::Float, "0.333"),
            Err(ParseError::InvalidNumberForType(
                Number::Count(0),
                Type::Float
            ))
        );

        assert_eq!(
            parse_value(Number::Count(1), Type::Float, "0.333"),
            Ok(Value::Float(0.333))
        );

        assert_eq!(
            parse_value(Number::Count(2), Type::Float, "0.333,0.667"),
            Ok(Value::Array(Array::Float(vec![Some(0.333), Some(0.667)])))
        );

        assert_eq!(
            parse_value(Number::Count(2), Type::Float, "0.333,."),
            Ok(Value::Array(Array::Float(vec![Some(0.333), None])))
        );
    }

    #[test]
    fn test_parse_value_with_character() {
        assert_eq!(
            parse_value(Number::Count(0), Type::Character, "n"),
            Err(ParseError::InvalidNumberForType(
                Number::Count(0),
                Type::Character
            ))
        );

        assert_eq!(
            parse_value(Number::Count(1), Type::Character, "n"),
            Ok(Value::Character('n'))
        );

        assert_eq!(
            parse_value(Number::Count(2), Type::Character, "n,d,l,s"),
            Ok(Value::Array(Array::Character(vec![
                Some('n'),
                Some('d'),
                Some('l'),
                Some('s')
            ])))
        );

        assert_eq!(
            parse_value(Number::Count(2), Type::Character, "n,d,l,."),
            Ok(Value::Array(Array::Character(vec![
                Some('n'),
                Some('d'),
                Some('l'),
                None
            ])))
        );
    }

    #[test]
    fn test_parse_value_with_string() {
        assert_eq!(
            parse_value(Number::Count(0), Type::String, "noodles"),
            Err(ParseError::InvalidNumberForType(
                Number::Count(0),
                Type::String
            ))
        );

        assert_eq!(
            parse_value(Number::Count(1), Type::String, "noodles"),
            Ok(Value::String(String::from("noodles")))
        );
        assert_eq!(
            parse_value(Number::Count(1), Type::String, "8%25"),
            Ok(Value::String(String::from("8%")))
        );

        assert_eq!(
            parse_value(Number::Count(2), Type::String, "noodles,vcf"),
            Ok(Value::Array(Array::String(vec![
                Some(String::from("noodles")),
                Some(String::from("vcf"))
            ])))
        );
        assert_eq!(
            parse_value(Number::Count(2), Type::String, "8%25,13%25"),
            Ok(Value::Array(Array::String(vec![
                Some(String::from("8%")),
                Some(String::from("13%")),
            ])))
        );
        assert_eq!(
            parse_value(Number::Count(2), Type::String, "noodles,."),
            Ok(Value::Array(Array::String(vec![
                Some(String::from("noodles")),
                None,
            ])))
        );
    }
}
