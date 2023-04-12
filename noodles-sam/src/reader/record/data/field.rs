mod tag;
mod ty;
mod value;

pub(crate) use self::value::parse_value;

use std::{error, fmt};

use self::{tag::parse_tag, ty::parse_type};
use crate::record::data::field::{Tag, Value};

/// An error when a raw BAM record data field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// Expected delimiter.
    ExpectedDelimiter,
    /// The tag is invalid.
    InvalidTag(tag::ParseError),
    /// The type is invalid.
    InvalidType(Tag, ty::ParseError),
    /// The value is invalid.
    InvalidValue(Tag, value::ParseError),
}

impl ParseError {
    /// Returns the tag of the field that caused the failure.
    pub fn tag(&self) -> Option<Tag> {
        match self {
            Self::InvalidType(tag, _) => Some(*tag),
            Self::InvalidValue(tag, _) => Some(*tag),
            _ => None,
        }
    }
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            ParseError::InvalidTag(e) => Some(e),
            ParseError::InvalidType(_, e) => Some(e),
            ParseError::InvalidValue(_, e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::UnexpectedEof => write!(f, "unexpected EOF"),
            ParseError::ExpectedDelimiter => write!(f, "expected delimiter"),
            ParseError::InvalidTag(_) => write!(f, "invalid tag"),
            ParseError::InvalidType(..) => write!(f, "invalid type"),
            ParseError::InvalidValue(..) => write!(f, "invalid value"),
        }
    }
}

pub(super) fn parse_field(src: &mut &[u8]) -> Result<(Tag, Value), ParseError> {
    use crate::reader::record::next_field;

    let mut buf = next_field(src);

    let tag = parse_tag(&mut buf).map_err(ParseError::InvalidTag)?;

    consume_delimiter(&mut buf)?;
    let ty = parse_type(&mut buf).map_err(|e| ParseError::InvalidType(tag, e))?;

    consume_delimiter(&mut buf)?;
    let value = parse_value(&mut buf, ty).map_err(|e| ParseError::InvalidValue(tag, e))?;

    Ok((tag, value))
}

fn consume_delimiter(src: &mut &[u8]) -> Result<(), ParseError> {
    const DELIMITER: u8 = b':';

    let (n, rest) = src.split_first().ok_or(ParseError::UnexpectedEof)?;

    *src = rest;

    if *n == DELIMITER {
        Ok(())
    } else {
        Err(ParseError::ExpectedDelimiter)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_field() {
        let mut src = &b"NH:i:1\tCO:Z:ndls"[..];

        assert_eq!(
            parse_field(&mut src),
            Ok((Tag::AlignmentHitCount, Value::from(1)))
        );

        assert_eq!(
            parse_field(&mut src),
            Ok((Tag::Comment, Value::String(String::from("ndls"))))
        );

        assert!(src.is_empty());
    }
}
