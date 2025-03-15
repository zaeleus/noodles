mod array;

use std::{borrow::Cow, io};

use self::array::Array;
use super::percent_decode;

/// A raw GFF record attributes field value.
#[derive(Debug, Eq, PartialEq)]
pub enum Value<'a> {
    /// A string.
    String(Cow<'a, str>),
    /// An array.
    Array(Array<'a>),
}

impl AsRef<str> for Value<'_> {
    fn as_ref(&self) -> &str {
        match self {
            Value::String(s) => s,
            Value::Array(array) => array.as_ref(),
        }
    }
}

impl<'a> From<Value<'a>> for crate::feature::record::attributes::field::Value<'a> {
    fn from(value: Value<'a>) -> Self {
        match value {
            Value::String(s) => Self::String(s),
            Value::Array(array) => Self::Array(Box::new(array)),
        }
    }
}

pub(super) fn parse_value(s: &str) -> io::Result<Value<'_>> {
    if is_array(s) {
        Ok(Value::Array(Array::new(s)))
    } else {
        percent_decode(s)
            .map(Value::String)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}

fn is_array(s: &str) -> bool {
    const SEPARATOR: char = ',';
    s.contains(SEPARATOR)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_value() -> io::Result<()> {
        assert_eq!(parse_value("ndls")?, Value::String(Cow::from("ndls")));
        assert_eq!(parse_value("8,13")?, Value::Array(Array::new("8,13")));
        Ok(())
    }

    #[test]
    fn test_is_array() {
        assert!(is_array("8,13"));
        assert!(!is_array("ndls"));
    }
}
