pub mod ty;

use std::{error, fmt, mem, str};

pub use self::ty::read_type;
use crate::lazy::record::{
    value::{Array, Float, Int16, Int32, Int8, Type},
    Value,
};

pub fn read_value<'a>(src: &mut &'a [u8]) -> Result<Option<Value<'a>>, DecodeError> {
    let ty = read_type(src).map_err(DecodeError::InvalidType)?;

    match ty {
        None => Ok(None),
        Some(Type::Int8(0)) => Ok(Some(Value::Int8(None))),
        Some(Type::Int8(1)) => read_i8(src),
        Some(Type::Int8(n)) => read_i8s(src, n),
        Some(Type::Int16(0)) => Ok(Some(Value::Int16(None))),
        Some(Type::Int16(1)) => read_i16(src),
        Some(Type::Int16(n)) => read_i16s(src, n),
        Some(Type::Int32(0)) => Ok(Some(Value::Int32(None))),
        Some(Type::Int32(1)) => read_i32(src),
        Some(Type::Int32(n)) => read_i32s(src, n),
        Some(Type::Float(0)) => Ok(Some(Value::Float(None))),
        Some(Type::Float(1)) => read_f32(src),
        Some(Type::Float(n)) => read_f32s(src, n),
        Some(Type::String(0)) => Ok(Some(Value::String(None))),
        Some(Type::String(n)) => read_string(src, n),
    }
}

fn read_i8<'a>(src: &mut &'a [u8]) -> Result<Option<Value<'a>>, DecodeError> {
    if let Some((b, rest)) = src.split_first() {
        *src = rest;
        Ok(Some(Value::Int8(Some(Int8::from(*b as i8)))))
    } else {
        Err(DecodeError::UnexpectedEof)
    }
}

fn read_i8s<'a>(src: &mut &'a [u8], len: usize) -> Result<Option<Value<'a>>, DecodeError> {
    if src.len() < len {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(len);
    let values = buf.iter().map(|&b| b as i8).collect();
    *src = rest;

    Ok(Some(Value::Array(Array::Int8(values))))
}

fn read_i16<'a>(src: &mut &'a [u8]) -> Result<Option<Value<'a>>, DecodeError> {
    if src.len() < mem::size_of::<i16>() {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(mem::size_of::<i16>());

    // SAFETY: `buf` is 2 bytes.
    let n = i16::from_le_bytes(buf.try_into().unwrap());

    *src = rest;

    Ok(Some(Value::Int16(Some(Int16::from(n)))))
}

fn read_i16s<'a>(src: &mut &'a [u8], len: usize) -> Result<Option<Value<'a>>, DecodeError> {
    let len = mem::size_of::<i16>() * len;

    if src.len() < len {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(len);

    let values = buf
        .chunks_exact(mem::size_of::<i16>())
        .map(|chunk| {
            // SAFETY: `chunk` is 2 bytes.
            i16::from_le_bytes(chunk.try_into().unwrap())
        })
        .collect();

    *src = rest;

    Ok(Some(Value::Array(Array::Int16(values))))
}

fn read_i32<'a>(src: &mut &'a [u8]) -> Result<Option<Value<'a>>, DecodeError> {
    if src.len() < mem::size_of::<i32>() {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(mem::size_of::<i32>());

    // SAFETY: `buf` is 4 bytes.
    let n = i32::from_le_bytes(buf.try_into().unwrap());

    *src = rest;

    Ok(Some(Value::Int32(Some(Int32::from(n)))))
}

fn read_i32s<'a>(src: &mut &'a [u8], len: usize) -> Result<Option<Value<'a>>, DecodeError> {
    let len = mem::size_of::<i32>() * len;

    if src.len() < len {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(len);

    let values = buf
        .chunks_exact(mem::size_of::<i32>())
        .map(|chunk| {
            // SAFETY: `chunk` is 4 bytes.
            i32::from_le_bytes(chunk.try_into().unwrap())
        })
        .collect();

    *src = rest;

    Ok(Some(Value::Array(Array::Int32(values))))
}

fn read_f32<'a>(src: &mut &'a [u8]) -> Result<Option<Value<'a>>, DecodeError> {
    if src.len() < mem::size_of::<f32>() {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(mem::size_of::<f32>());

    // SAFETY: `buf` is 4 bytes.
    let n = f32::from_le_bytes(buf.try_into().unwrap());

    *src = rest;

    Ok(Some(Value::Float(Some(Float::from(n)))))
}

fn read_f32s<'a>(src: &mut &'a [u8], len: usize) -> Result<Option<Value<'a>>, DecodeError> {
    let len = mem::size_of::<f32>() * len;

    if src.len() < len {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(len);

    let values = buf
        .chunks_exact(mem::size_of::<f32>())
        .map(|chunk| {
            // SAFETY: `chunk` is 4 bytes.
            f32::from_le_bytes(chunk.try_into().unwrap())
        })
        .collect();

    *src = rest;

    Ok(Some(Value::Array(Array::Float(values))))
}

fn read_string<'a>(src: &mut &'a [u8], len: usize) -> Result<Option<Value<'a>>, DecodeError> {
    if src.len() < len {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(len);
    let s = str::from_utf8(buf).map_err(DecodeError::InvalidString)?;
    *src = rest;

    Ok(Some(Value::String(Some(s))))
}

#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    UnexpectedEof,
    InvalidType(ty::DecodeError),
    InvalidString(str::Utf8Error),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::UnexpectedEof => None,
            Self::InvalidType(e) => Some(e),
            Self::InvalidString(e) => Some(e),
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidType(_) => write!(f, "invalid type"),
            Self::InvalidString(_) => write!(f, "invalid string"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_value() {
        fn t(mut src: &[u8], expected: Option<Value<'_>>) {
            assert_eq!(read_value(&mut src), Ok(expected));
        }

        t(&[0x00], None);
        t(&[0x01], Some(Value::Int8(None)));
        t(&[0x11, 0x05], Some(Value::Int8(Some(Int8::Value(5)))));

        let src = &[0x31, 0x05, 0x08, 0x0d];
        t(src, Some(Value::Array(Array::Int8(vec![5, 8, 13]))));

        t(&[0x02], Some(Value::Int16(None)));

        let src = &[0x12, 0x79, 0x01];
        t(src, Some(Value::Int16(Some(Int16::Value(377)))));

        let src = &[0x32, 0x79, 0x01, 0x62, 0x02, 0xdb, 0x03];
        t(src, Some(Value::Array(Array::Int16(vec![377, 610, 987]))));

        t(&[0x03], Some(Value::Int32(None)));

        let src = &[0x13, 0x11, 0x25, 0x01, 0x00];
        t(src, Some(Value::Int32(Some(Int32::Value(75025)))));

        let src = &[
            0x33, 0x11, 0x25, 0x01, 0x00, 0x31, 0xda, 0x01, 0x00, 0x42, 0xff, 0x02, 0x00,
        ];
        t(
            src,
            Some(Value::Array(Array::Int32(vec![75025, 121393, 196418]))),
        );

        t(&[0x05], Some(Value::Float(None)));

        let src = &[0x15, 0x00, 0x00, 0x00, 0x00];
        t(src, Some(Value::Float(Some(Float::from(0.0)))));

        let src = &[0x25, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3f];
        t(src, Some(Value::Array(Array::Float(vec![0.0, 0.5]))));

        t(&[0x07], Some(Value::String(None)));
        t(&[0x17, b'n'], Some(Value::String(Some("n"))));

        let src = &[0x47, b'n', b'd', b'l', b's'];
        t(src, Some(Value::String(Some("ndls"))));
    }
}
