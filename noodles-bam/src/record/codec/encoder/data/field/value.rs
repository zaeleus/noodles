pub mod array;

use std::io;

use bytes::BufMut;
use noodles_sam::alignment::record::data::field::Value;

use self::array::put_array;

pub fn put_value<B>(dst: &mut B, value: &Value) -> io::Result<()>
where
    B: BufMut,
{
    match value {
        Value::Character(c) => dst.put_u8(*c),
        Value::Int8(n) => dst.put_i8(*n),
        Value::UInt8(n) => dst.put_u8(*n),
        Value::Int16(n) => dst.put_i16_le(*n),
        Value::UInt16(n) => dst.put_u16_le(*n),
        Value::Int32(n) => dst.put_i32_le(*n),
        Value::UInt32(n) => dst.put_u32_le(*n),
        Value::Float(n) => dst.put_f32_le(*n),
        Value::String(s) => put_string(dst, s),
        Value::Hex(s) => put_string(dst, s),
        Value::Array(array) => put_array(dst, array)?,
    }

    Ok(())
}

fn put_string<B>(dst: &mut B, buf: &[u8])
where
    B: BufMut,
{
    const NUL: u8 = 0x00;

    dst.put(buf);
    dst.put_u8(NUL);
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;

    use super::*;

    #[test]
    fn test_put_value() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_sam::alignment::record::data::field::value::Array;

        fn t(buf: &mut Vec<u8>, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            put_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Value::Character(b'n'), &[b'n'])?;
        t(&mut buf, &Value::Int8(1), &[0x01])?;
        t(&mut buf, &Value::UInt8(2), &[0x02])?;
        t(&mut buf, &Value::Int16(3), &[0x03, 0x00])?;
        t(&mut buf, &Value::UInt16(5), &[0x05, 0x00])?;
        t(&mut buf, &Value::Int32(8), &[0x08, 0x00, 0x00, 0x00])?;
        t(&mut buf, &Value::UInt32(13), &[0x0d, 0x00, 0x00, 0x00])?;
        t(&mut buf, &Value::Float(8.0), &[0x00, 0x00, 0x00, 0x41])?;

        t(
            &mut buf,
            &Value::String(b"ndls".as_bstr()),
            &[b'n', b'd', b'l', b's', 0x00],
        )?;

        t(
            &mut buf,
            &Value::Hex(b"CAFE".as_bstr()),
            &[b'C', b'A', b'F', b'E', 0x00],
        )?;

        use super::array::tests::T;

        t(
            &mut buf,
            &Value::Array(Array::UInt8(Box::new(T::new(&[0])))),
            &[
                b'C', // subtype = UInt8
                0x01, 0x00, 0x00, 0x00, // count = 2
                0x00, // values[0] = 0
            ],
        )?;

        Ok(())
    }
}
