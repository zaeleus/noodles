pub mod subtype;

use std::{error, fmt, num};

use noodles_sam::alignment::{
    record::data::field::value::array::Subtype,
    record_buf::data::field::{Value, value::Array},
};

use self::subtype::read_subtype;

// An error when a raw BAM record data field array value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The subtype is invalid.
    InvalidSubtype(subtype::DecodeError),
    /// The length is invalid.
    InvalidLength(num::TryFromIntError),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::UnexpectedEof => None,
            Self::InvalidSubtype(e) => Some(e),
            Self::InvalidLength(e) => Some(e),
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidSubtype(_) => write!(f, "invalid subtype"),
            Self::InvalidLength(_) => write!(f, "invalid length"),
        }
    }
}

pub(super) fn get_array(src: &mut &[u8]) -> Result<Value, DecodeError> {
    let subtype = read_subtype(src).map_err(DecodeError::InvalidSubtype)?;

    let len =
        read_u32_le(src).and_then(|n| usize::try_from(n).map_err(DecodeError::InvalidLength))?;

    let array = match subtype {
        Subtype::Int8 => Array::Int8((0..len).map(|_| read_i8(src)).collect::<Result<_, _>>()?),
        Subtype::UInt8 => Array::UInt8((0..len).map(|_| read_u8(src)).collect::<Result<_, _>>()?),
        Subtype::Int16 => Array::Int16(
            (0..len)
                .map(|_| read_i16_le(src))
                .collect::<Result<_, _>>()?,
        ),

        Subtype::UInt16 => Array::UInt16(
            (0..len)
                .map(|_| read_u16_le(src))
                .collect::<Result<_, _>>()?,
        ),

        Subtype::Int32 => Array::Int32(
            (0..len)
                .map(|_| read_i32_le(src))
                .collect::<Result<_, _>>()?,
        ),

        Subtype::UInt32 => Array::UInt32(
            (0..len)
                .map(|_| read_u32_le(src))
                .collect::<Result<_, _>>()?,
        ),

        Subtype::Float => Array::Float(
            (0..len)
                .map(|_| read_f32_le(src))
                .collect::<Result<_, _>>()?,
        ),
    };

    Ok(Value::Array(array))
}

fn read_i8(src: &mut &[u8]) -> Result<i8, DecodeError> {
    let (n, rest) = src.split_first().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(*n as i8)
}

fn read_u8(src: &mut &[u8]) -> Result<u8, DecodeError> {
    let (n, rest) = src.split_first().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(*n)
}

fn read_i16_le(src: &mut &[u8]) -> Result<i16, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(i16::from_le_bytes(*buf))
}

fn read_u16_le(src: &mut &[u8]) -> Result<u16, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(u16::from_le_bytes(*buf))
}

fn read_i32_le(src: &mut &[u8]) -> Result<i32, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(i32::from_le_bytes(*buf))
}

fn read_u32_le(src: &mut &[u8]) -> Result<u32, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(u32::from_le_bytes(*buf))
}

fn read_f32_le(src: &mut &[u8]) -> Result<f32, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(f32::from_le_bytes(*buf))
}
