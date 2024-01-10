pub mod subtype;

use std::{error, fmt, num};

use bytes::Buf;
use noodles_sam::{
    alignment::{
        record::data::field::value::array::Subtype, record_buf::data::field::value::Array,
    },
    record::data::field::Value,
};

use self::subtype::get_subtype;

// An error when a raw BAM record data field array value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// The subtype is invalid.
    InvalidSubtype(subtype::DecodeError),
    /// The length is invalid.
    InvalidLength(num::TryFromIntError),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidSubtype(e) => Some(e),
            Self::InvalidLength(e) => Some(e),
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidSubtype(_) => write!(f, "invalid subtype"),
            Self::InvalidLength(_) => write!(f, "invalid length"),
        }
    }
}

pub(super) fn get_array<B>(src: &mut B) -> Result<Value, DecodeError>
where
    B: Buf,
{
    let subtype = get_subtype(src).map_err(DecodeError::InvalidSubtype)?;
    let len = usize::try_from(src.get_u32_le()).map_err(DecodeError::InvalidLength)?;

    match subtype {
        Subtype::Int8 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_i8());
            }

            Ok(Value::Array(Array::Int8(buf)))
        }
        Subtype::UInt8 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_u8());
            }

            Ok(Value::Array(Array::UInt8(buf)))
        }
        Subtype::Int16 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_i16_le());
            }

            Ok(Value::Array(Array::Int16(buf)))
        }
        Subtype::UInt16 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_u16_le());
            }

            Ok(Value::Array(Array::UInt16(buf)))
        }
        Subtype::Int32 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_i32_le());
            }

            Ok(Value::Array(Array::Int32(buf)))
        }
        Subtype::UInt32 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_u32_le());
            }

            Ok(Value::Array(Array::UInt32(buf)))
        }
        Subtype::Float => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_f32_le());
            }

            Ok(Value::Array(Array::Float(buf)))
        }
    }
}
