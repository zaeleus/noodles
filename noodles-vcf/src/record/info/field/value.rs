//! VCF record info field value.

mod array;

use std::{fmt, str};

pub use self::array::Array;

/// A VCF record info field value.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    /// An 32-bit integer.
    Integer(i32),
    /// A single-precision floating-point.
    Float(f32),
    /// A boolean.
    Flag,
    /// A character.
    Character(char),
    /// A string.
    String(String),
    /// An array.
    Array(Array),
}

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Integer(n) => write!(f, "{n}"),
            Self::Float(n) => write!(f, "{n}"),
            Self::Flag => Ok(()),
            Self::Character(c) => write!(f, "{c}"),
            Self::String(s) => write!(f, "{s}"),
            Self::Array(array) => write!(f, "{array}"),
        }
    }
}

impl From<i32> for Value {
    fn from(n: i32) -> Self {
        Self::Integer(n)
    }
}

impl From<f32> for Value {
    fn from(n: f32) -> Self {
        Self::Float(n)
    }
}

impl From<char> for Value {
    fn from(c: char) -> Self {
        Self::Character(c)
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

impl From<Vec<Option<i32>>> for Value {
    fn from(values: Vec<Option<i32>>) -> Self {
        Self::Array(Array::Integer(values))
    }
}

impl From<Vec<Option<f32>>> for Value {
    fn from(values: Vec<Option<f32>>) -> Self {
        Self::Array(Array::Float(values))
    }
}

impl From<Vec<Option<char>>> for Value {
    fn from(values: Vec<Option<char>>) -> Self {
        Self::Array(Array::Character(values))
    }
}

impl From<Vec<Option<String>>> for Value {
    fn from(values: Vec<Option<String>>) -> Self {
        Self::Array(Array::String(values))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let value = Value::from(2);
        assert_eq!(value.to_string(), "2");

        let value = Value::from(0.333);
        assert_eq!(value.to_string(), "0.333");

        assert_eq!(Value::Flag.to_string(), "");

        let value = Value::from('n');
        assert_eq!(value.to_string(), "n");

        let value = Value::from("noodles");
        assert_eq!(value.to_string(), "noodles");

        let value = Value::from(vec![Some(2)]);
        assert_eq!(value.to_string(), "2");

        let value = Value::from(vec![Some(2), Some(5)]);
        assert_eq!(value.to_string(), "2,5");

        let value = Value::from(vec![Some(2), None]);
        assert_eq!(value.to_string(), "2,.");

        let value = Value::from(vec![Some(0.333)]);
        assert_eq!(value.to_string(), "0.333");

        let value = Value::from(vec![Some(0.333), Some(0.667)]);
        assert_eq!(value.to_string(), "0.333,0.667");

        let value = Value::from(vec![Some(0.333), None]);
        assert_eq!(value.to_string(), "0.333,.");

        let value = Value::from(vec![Some('n')]);
        assert_eq!(value.to_string(), "n");

        let value = Value::from(vec![Some('n'), Some('d'), Some('l'), Some('s')]);
        assert_eq!(value.to_string(), "n,d,l,s");

        let value = Value::from(vec![Some('n'), Some('d'), Some('l'), None]);
        assert_eq!(value.to_string(), "n,d,l,.");

        let value = Value::from(vec![Some(String::from("noodles"))]);
        assert_eq!(value.to_string(), "noodles");

        let value = Value::from(vec![
            Some(String::from("noodles")),
            Some(String::from("vcf")),
        ]);
        assert_eq!(value.to_string(), "noodles,vcf");

        let value = Value::from(vec![Some(String::from("noodles")), None]);
        assert_eq!(value.to_string(), "noodles,.");
    }
}
