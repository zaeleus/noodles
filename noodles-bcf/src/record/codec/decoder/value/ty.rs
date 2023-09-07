use std::{error, fmt, num};

use super::read_value;
use crate::lazy::record::value::Type;

pub fn read_type(src: &mut &[u8]) -> Result<Option<Type>, DecodeError> {
    let encoding = get_u8(src)?;
    let mut len = usize::from(encoding >> 4);

    if len == 0x0f {
        let value = read_value(src).map_err(|e| DecodeError::InvalidValue(Box::new(e)))?;

        len = match value.and_then(|v| v.as_int()) {
            Some(n) => usize::try_from(n).map_err(DecodeError::InvalidLength)?,
            None => return Err(DecodeError::InvalidLengthValue),
        };
    }

    match encoding & 0x0f {
        0 => Ok(None),
        1 => Ok(Some(Type::Int8(len))),
        2 => Ok(Some(Type::Int16(len))),
        3 => Ok(Some(Type::Int32(len))),
        5 => Ok(Some(Type::Float(len))),
        7 => Ok(Some(Type::String(len))),
        ty => Err(DecodeError::InvalidType(ty)),
    }
}

fn get_u8(src: &mut &[u8]) -> Result<u8, DecodeError> {
    if let Some((b, rest)) = src.split_first() {
        *src = rest;
        Ok(*b)
    } else {
        Err(DecodeError::UnexpectedEof)
    }
}

#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    UnexpectedEof,
    InvalidValue(Box<super::DecodeError>),
    InvalidLength(num::TryFromIntError),
    InvalidLengthValue,
    InvalidType(u8),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidLength(e) => Some(e),
            Self::InvalidValue(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidValue(_) => write!(f, "invalid value"),
            Self::InvalidLength(_) => write!(f, "invalid length"),
            Self::InvalidLengthValue => write!(f, "invalid length value"),
            Self::InvalidType(actual) => write!(
                f,
                "invalid type: expected {{0, 1, 2, 3, 5, 7}}, got {actual}"
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_type() {
        fn t(mut src: &[u8], expected: Option<Type>) {
            assert_eq!(read_type(&mut src), Ok(expected));
        }

        t(&[0x00], None);
        t(&[0x10], None);
        t(&[0x11], Some(Type::Int8(1)));
        t(&[0x12], Some(Type::Int16(1)));
        t(&[0x13], Some(Type::Int32(1)));
        t(&[0x15], Some(Type::Float(1)));
        t(&[0x17], Some(Type::String(1)));

        t(
            &[
                0xf1, // (len >= 15, Type::Int8)
                0x11, // Some(Type::Int8(1))
                0x15, // Some(Value::Int8(21))
            ],
            Some(Type::Int8(21)),
        );

        let mut src = &[][..];
        assert_eq!(read_type(&mut src), Err(DecodeError::UnexpectedEof));

        let mut src = &[0x14][..];
        assert_eq!(read_type(&mut src), Err(DecodeError::InvalidType(4)));
    }
}
