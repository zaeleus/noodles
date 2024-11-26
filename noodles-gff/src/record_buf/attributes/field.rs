//! GFF record attributes field.

pub mod tag;
pub mod value;

use std::{borrow::Cow, error, fmt, str};

use percent_encoding::percent_decode_str;

pub use self::{tag::Tag, value::Value};

pub(super) const SEPARATOR: char = '=';

/// An error returned when a raw attributes field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
    /// A key is invalid.
    InvalidKey(str::Utf8Error),
    /// A value is invalid.
    InvalidValue(Tag, value::ParseError),
}

impl ParseError {
    /// Returns the key of the field that caused the failure.
    pub fn key(&self) -> Option<&Tag> {
        match self {
            Self::InvalidValue(key, _) => Some(key),
            _ => None,
        }
    }
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKey(e) => Some(e),
            Self::InvalidValue(_, e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => write!(f, "invalid input"),
            Self::InvalidKey(_) => write!(f, "invalid key"),
            Self::InvalidValue(..) => write!(f, "invalid value"),
        }
    }
}

pub(super) fn parse_field(s: &str) -> Result<(Tag, Value), ParseError> {
    let (raw_key, raw_value) = s.split_once(SEPARATOR).ok_or(ParseError::Invalid)?;

    let key: Tag = percent_decode(raw_key)
        .map(|k| k.into())
        .map_err(ParseError::InvalidKey)?;

    let value = raw_value
        .parse()
        .map_err(|e| ParseError::InvalidValue(key.clone(), e))?;

    Ok((key, value))
}

pub(super) fn percent_decode(s: &str) -> Result<Cow<'_, str>, str::Utf8Error> {
    percent_decode_str(s).decode_utf8()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_field() {
        assert_eq!(
            parse_field("gene_id=gene0"),
            Ok((String::from("gene_id"), Value::from("gene0")))
        );

        assert_eq!(
            parse_field("%25s=13%2C21"),
            Ok((String::from("%s"), Value::from("13,21")))
        );

        assert_eq!(parse_field(""), Err(ParseError::Invalid));
    }
}
