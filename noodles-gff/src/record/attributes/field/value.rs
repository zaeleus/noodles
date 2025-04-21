mod array;

use std::borrow::Cow;

use bstr::BStr;

use self::array::Array;
use super::percent_decode;

/// A raw GFF record attributes field value.
#[derive(Debug, Eq, PartialEq)]
pub enum Value<'a> {
    /// A string.
    String(Cow<'a, BStr>),
    /// An array.
    Array(Array<'a>),
}

impl AsRef<BStr> for Value<'_> {
    fn as_ref(&self) -> &BStr {
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

pub(crate) fn parse_value(src: &[u8]) -> Value<'_> {
    if is_array(src) {
        Value::Array(Array::new(src))
    } else {
        Value::String(percent_decode(src))
    }
}

fn is_array(src: &[u8]) -> bool {
    const SEPARATOR: u8 = b',';
    src.contains(&SEPARATOR)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_value() {
        assert_eq!(
            parse_value(b"ndls"),
            Value::String(Cow::from(BStr::new("ndls")))
        );
        assert_eq!(parse_value(b"8,13"), Value::Array(Array::new(b"8,13")));
    }

    #[test]
    fn test_is_array() {
        assert!(is_array(b"8,13"));
        assert!(!is_array(b"ndls"));
    }
}
