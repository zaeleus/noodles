//! BAM record data field component readers.

mod tag;
mod ty;
mod value;

pub use self::value::get_value;

use std::{error, fmt};

use bytes::Buf;
use noodles_sam::record::data::field::{Tag, Value};

use self::{tag::get_tag, ty::get_type};

/// An error when a raw BAM record data field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
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
            Self::InvalidTag(_) => None,
            Self::InvalidType(tag, _) => Some(*tag),
            Self::InvalidValue(tag, _) => Some(*tag),
        }
    }
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            ParseError::InvalidTag(e) => Some(e),
            ParseError::InvalidType(_, e) => Some(e),
            ParseError::InvalidValue(_, e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::InvalidTag(_) => write!(f, "invalid tag"),
            ParseError::InvalidType(..) => write!(f, "invalid type"),
            ParseError::InvalidValue(..) => write!(f, "invalid value"),
        }
    }
}

pub(crate) fn get_field<B>(src: &mut B) -> Result<(Tag, Value), ParseError>
where
    B: Buf,
{
    let tag = get_tag(src).map_err(ParseError::InvalidTag)?;

    let ty = get_type(src).map_err(|e| ParseError::InvalidType(tag, e))?;
    let value = get_value(src, ty).map_err(|e| ParseError::InvalidValue(tag, e))?;

    Ok((tag, value))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_field() {
        use noodles_sam::record::data::field::{tag, Value};

        let data = [b'N', b'H', b'C', 0x01];
        let mut reader = &data[..];
        assert_eq!(
            get_field(&mut reader),
            Ok((tag::ALIGNMENT_HIT_COUNT, Value::from(1)))
        );

        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            get_field(&mut reader),
            Err(ParseError::InvalidTag(_))
        ));

        let data = [b'N', b'H', b'z'];
        let mut reader = &data[..];
        assert!(matches!(
            get_field(&mut reader),
            Err(ParseError::InvalidType(tag::ALIGNMENT_HIT_COUNT, _))
        ));

        let data = [b'N', b'H', b'C'];
        let mut reader = &data[..];
        assert!(matches!(
            get_field(&mut reader),
            Err(ParseError::InvalidValue(tag::ALIGNMENT_HIT_COUNT, _))
        ));
    }
}
