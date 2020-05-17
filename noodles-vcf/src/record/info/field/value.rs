use std::{error, fmt, num};

use crate::header::{info::Type, Number};

use super::Key;

const DELIMITER: char = ',';

#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    Integer(i32),
    Float(f32),
    Flag,
    Character(char),
    String(String),
    IntegerArray(Vec<i32>),
    FloatArray(Vec<f32>),
    CharacterArray(Vec<char>),
    StringArray(Vec<String>),
}

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::Integer(n) => write!(f, "{}", n),
            Self::Float(n) => write!(f, "{}", n),
            Self::Flag => Ok(()),
            Self::Character(c) => write!(f, "{}", c),
            Self::String(s) => write!(f, "{}", s),
            Self::IntegerArray(values) => {
                for (i, value) in values.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{}", DELIMITER)?;
                    }

                    write!(f, "{}", value)?;
                }

                Ok(())
            }
            Self::FloatArray(values) => {
                for (i, value) in values.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{}", DELIMITER)?;
                    }

                    write!(f, "{}", value)?;
                }

                Ok(())
            }
            Self::CharacterArray(values) => {
                for (i, value) in values.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{}", DELIMITER)?;
                    }

                    write!(f, "{}", value)?;
                }

                Ok(())
            }
            Self::StringArray(values) => {
                for (i, value) in values.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{}", DELIMITER)?;
                    }

                    write!(f, "{}", value)?;
                }

                Ok(())
            }
        }
    }
}

#[derive(Debug)]
pub enum ParseError {
    InvalidNumberForType(Number, Type),
    InvalidInteger(num::ParseIntError),
    InvalidFloat(num::ParseFloatError),
    InvalidFlag(String),
    InvalidCharacter(String),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("invalid info value: ")?;

        match self {
            ParseError::InvalidNumberForType(number, ty) => {
                write!(f, "invalid number {:?} for type {:?}", number, ty)
            }
            ParseError::InvalidInteger(e) => write!(f, "invalid integer: {}", e),
            ParseError::InvalidFloat(e) => write!(f, "invalid float: {}", e),
            ParseError::InvalidFlag(e) => write!(f, "invalid flag: {}", e),
            ParseError::InvalidCharacter(e) => write!(f, "invalid character: {}", e),
        }
    }
}

impl Value {
    pub fn from_str_key(s: &str, key: &Key) -> Result<Self, ParseError> {
        match key.ty() {
            Type::Integer => match key.number() {
                Number::Count(0) => Err(ParseError::InvalidNumberForType(key.number(), key.ty())),
                Number::Count(1) => parse_i32(s),
                _ => parse_i32_array(s),
            },
            Type::Float => match key.number() {
                Number::Count(0) => Err(ParseError::InvalidNumberForType(key.number(), key.ty())),
                Number::Count(1) => parse_f32(s),
                _ => parse_f32_array(s),
            },
            Type::Flag => match key.number() {
                Number::Count(0) => parse_flag(s),
                _ => Err(ParseError::InvalidNumberForType(key.number(), key.ty())),
            },
            Type::Character => match key.number() {
                Number::Count(0) => Err(ParseError::InvalidNumberForType(key.number(), key.ty())),
                Number::Count(1) => parse_char(s),
                _ => parse_char_array(s),
            },
            Type::String => match key.number() {
                Number::Count(0) => Err(ParseError::InvalidNumberForType(key.number(), key.ty())),
                Number::Count(1) => Ok(parse_string(s)),
                _ => Ok(parse_string_array(s)),
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
        .map(|s| s.parse().map_err(ParseError::InvalidInteger))
        .collect::<Result<_, _>>()
        .map(Value::IntegerArray)
}

fn parse_f32(s: &str) -> Result<Value, ParseError> {
    s.parse()
        .map(Value::Float)
        .map_err(ParseError::InvalidFloat)
}

fn parse_f32_array(s: &str) -> Result<Value, ParseError> {
    s.split(DELIMITER)
        .map(|s| s.parse().map_err(ParseError::InvalidFloat))
        .collect::<Result<_, _>>()
        .map(Value::FloatArray)
}

fn parse_flag(s: &str) -> Result<Value, ParseError> {
    if s.is_empty() {
        Ok(Value::Flag)
    } else {
        Err(ParseError::InvalidFlag(s.into()))
    }
}

fn parse_raw_char(s: &str) -> Result<char, ParseError> {
    let mut chars = s.chars();

    if let Some(c) = chars.next() {
        if chars.next().is_none() {
            return Ok(c);
        }
    }

    Err(ParseError::InvalidCharacter(s.into()))
}

fn parse_char(s: &str) -> Result<Value, ParseError> {
    parse_raw_char(s).map(Value::Character)
}

fn parse_char_array(s: &str) -> Result<Value, ParseError> {
    s.split(DELIMITER)
        .map(parse_raw_char)
        .collect::<Result<_, _>>()
        .map(Value::CharacterArray)
}

fn parse_string(s: &str) -> Value {
    Value::String(s.into())
}

fn parse_string_array(s: &str) -> Value {
    let values = s.split(DELIMITER).map(|t| t.into()).collect();
    Value::StringArray(values)
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

        assert_eq!(Value::Flag.to_string(), "");

        let value = Value::Character('n');
        assert_eq!(value.to_string(), "n");

        let value = Value::String(String::from("noodles"));
        assert_eq!(value.to_string(), "noodles");

        let value = Value::IntegerArray(vec![2]);
        assert_eq!(value.to_string(), "2");

        let value = Value::IntegerArray(vec![2, 5]);
        assert_eq!(value.to_string(), "2,5");

        let value = Value::FloatArray(vec![0.333]);
        assert_eq!(value.to_string(), "0.333");

        let value = Value::FloatArray(vec![0.333, 0.667]);
        assert_eq!(value.to_string(), "0.333,0.667");

        let value = Value::CharacterArray(vec!['n']);
        assert_eq!(value.to_string(), "n");

        let value = Value::CharacterArray(vec!['n', 'd', 'l', 's']);
        assert_eq!(value.to_string(), "n,d,l,s");

        let value = Value::StringArray(vec![String::from("noodles")]);
        assert_eq!(value.to_string(), "noodles");

        let value = Value::StringArray(vec![String::from("noodles"), String::from("vcf")]);
        assert_eq!(value.to_string(), "noodles,vcf");
    }
}
