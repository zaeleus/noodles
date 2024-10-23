mod tag;
mod ty;
mod value;

use std::{error, fmt};

use noodles_sam::alignment::{record::data::field::Tag, record_buf::data::field::Value};

pub use self::value::read_value;
use self::{tag::read_tag, ty::read_type};

/// An error when a raw BAM record data field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// The tag is invalid.
    InvalidTag(tag::DecodeError),
    /// The type is invalid.
    InvalidType(Tag, ty::DecodeError),
    /// The value is invalid.
    InvalidValue(Tag, value::DecodeError),
}

impl DecodeError {
    /// Returns the tag of the field that caused the failure.
    pub fn tag(&self) -> Option<Tag> {
        match self {
            Self::InvalidTag(_) => None,
            Self::InvalidType(tag, _) => Some(*tag),
            Self::InvalidValue(tag, _) => Some(*tag),
        }
    }
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            DecodeError::InvalidTag(e) => Some(e),
            DecodeError::InvalidType(_, e) => Some(e),
            DecodeError::InvalidValue(_, e) => Some(e),
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DecodeError::InvalidTag(_) => write!(f, "invalid tag"),
            DecodeError::InvalidType(..) => write!(f, "invalid type"),
            DecodeError::InvalidValue(..) => write!(f, "invalid value"),
        }
    }
}

pub(crate) fn read_field(src: &mut &[u8]) -> Result<(Tag, Value), DecodeError> {
    let tag = read_tag(src).map_err(DecodeError::InvalidTag)?;

    let ty = read_type(src).map_err(|e| DecodeError::InvalidType(tag, e))?;
    let value = read_value(src, ty).map_err(|e| DecodeError::InvalidValue(tag, e))?;

    Ok((tag, value))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_field() {
        let data = [b'N', b'H', b'C', 0x01];
        let mut reader = &data[..];
        assert_eq!(
            read_field(&mut reader),
            Ok((Tag::ALIGNMENT_HIT_COUNT, Value::from(1)))
        );

        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_field(&mut reader),
            Err(DecodeError::InvalidTag(_))
        ));

        let data = [b'N', b'H', b'z'];
        let mut reader = &data[..];
        assert!(matches!(
            read_field(&mut reader),
            Err(DecodeError::InvalidType(Tag::ALIGNMENT_HIT_COUNT, _))
        ));

        let data = [b'N', b'H', b'C'];
        let mut reader = &data[..];
        assert!(matches!(
            read_field(&mut reader),
            Err(DecodeError::InvalidValue(Tag::ALIGNMENT_HIT_COUNT, _))
        ));
    }
}
