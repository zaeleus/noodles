//! BAM record data field component readers.

mod tag;
mod value;

pub use self::value::get_value;

use std::{error, fmt};

use bytes::Buf;
use noodles_sam::record::data::field::{Tag, Value};

/// An error when a raw BAM record data field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The tag is invalid.
    InvalidTag(tag::ParseError),
    /// The type is invalid.
    InvalidType(value::ty::ParseError),
    /// The value is invalid.
    InvalidValue(value::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            ParseError::InvalidTag(e) => Some(e),
            ParseError::InvalidType(e) => Some(e),
            ParseError::InvalidValue(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::InvalidTag(_) => write!(f, "invalid tag"),
            ParseError::InvalidType(_) => write!(f, "invalid type"),
            ParseError::InvalidValue(_) => write!(f, "invalid value"),
        }
    }
}

pub(crate) fn get_field<B>(src: &mut B) -> Result<(Tag, Value), ParseError>
where
    B: Buf,
{
    use self::tag::get_tag;

    let tag = get_tag(src).map_err(ParseError::InvalidTag)?;

    let ty = value::get_type(src).map_err(ParseError::InvalidType)?;
    let value = get_value(src, ty).map_err(ParseError::InvalidValue)?;

    Ok((tag, value))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_field() {
        use noodles_sam::record::data::field::{Tag, Value};

        let data = [b'N', b'H', b'C', 0x01];
        let mut reader = &data[..];
        assert_eq!(
            get_field(&mut reader),
            Ok((Tag::AlignmentHitCount, Value::from(1)))
        );

        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            get_field(&mut reader),
            Err(ParseError::InvalidTag(_))
        ));
    }
}
