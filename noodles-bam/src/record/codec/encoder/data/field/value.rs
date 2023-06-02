pub mod array;

use std::{ffi::CString, io};

use bytes::BufMut;
use noodles_sam::record::data::field::Value;

use self::array::put_array;

pub fn put_value<B>(dst: &mut B, value: &Value) -> io::Result<()>
where
    B: BufMut,
{
    match value {
        Value::Character(c) => dst.put_u8(u8::from(*c)),
        Value::Int8(n) => dst.put_i8(*n),
        Value::UInt8(n) => dst.put_u8(*n),
        Value::Int16(n) => dst.put_i16_le(*n),
        Value::UInt16(n) => dst.put_u16_le(*n),
        Value::Int32(n) => dst.put_i32_le(*n),
        Value::UInt32(n) => dst.put_u32_le(*n),
        Value::Float(n) => dst.put_f32_le(*n),
        Value::String(s) => put_string(dst, s)?,
        Value::Hex(s) => put_string(dst, s.as_ref())?,
        Value::Array(array) => put_array(dst, array)?,
    }

    Ok(())
}

fn put_string<B>(dst: &mut B, s: &str) -> io::Result<()>
where
    B: BufMut,
{
    let c_str = CString::new(s).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.put(c_str.as_bytes_with_nul());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_value() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_sam::record::data::field::value::{Array, Character};

        fn t(buf: &mut Vec<u8>, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            put_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(
            &mut buf,
            &Value::Character(Character::try_from('n')?),
            &[b'n'],
        )?;
        t(&mut buf, &Value::Int8(1), &[0x01])?;
        t(&mut buf, &Value::UInt8(2), &[0x02])?;
        t(&mut buf, &Value::Int16(3), &[0x03, 0x00])?;
        t(&mut buf, &Value::UInt16(5), &[0x05, 0x00])?;
        t(&mut buf, &Value::Int32(8), &[0x08, 0x00, 0x00, 0x00])?;
        t(&mut buf, &Value::UInt32(13), &[0x0d, 0x00, 0x00, 0x00])?;
        t(&mut buf, &Value::Float(8.0), &[0x00, 0x00, 0x00, 0x41])?;

        t(
            &mut buf,
            &Value::String(String::from("ndls")),
            &[b'n', b'd', b'l', b's', 0x00],
        )?;

        t(
            &mut buf,
            &Value::Hex("CAFE".parse()?),
            &[b'C', b'A', b'F', b'E', 0x00],
        )?;

        t(
            &mut buf,
            &Value::Array(Array::UInt8(vec![0])),
            &[
                b'C', // subtype = UInt8
                0x01, 0x00, 0x00, 0x00, // count = 2
                0x00, // values[0] = 0
            ],
        )?;

        Ok(())
    }
}
