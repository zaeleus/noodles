pub mod subtype;

use bytes::Buf;
use noodles_sam::record::data::field::{value::Subtype, Value};

use super::{get_subtype, ParseError};

pub(super) fn get_array_value<B>(src: &mut B) -> Result<Value, ParseError>
where
    B: Buf,
{
    let subtype = get_subtype(src).map_err(ParseError::InvalidArraySubtype)?;
    let len = usize::try_from(src.get_i32_le()).map_err(ParseError::InvalidArrayLength)?;

    match subtype {
        Subtype::Int8 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_i8());
            }

            Ok(Value::Int8Array(buf))
        }
        Subtype::UInt8 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_u8());
            }

            Ok(Value::UInt8Array(buf))
        }
        Subtype::Int16 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_i16_le());
            }

            Ok(Value::Int16Array(buf))
        }
        Subtype::UInt16 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_u16_le());
            }

            Ok(Value::UInt16Array(buf))
        }
        Subtype::Int32 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_i32_le());
            }

            Ok(Value::Int32Array(buf))
        }
        Subtype::UInt32 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_u32_le());
            }

            Ok(Value::UInt32Array(buf))
        }
        Subtype::Float => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_f32_le());
            }

            Ok(Value::FloatArray(buf))
        }
    }
}
