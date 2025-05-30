pub mod array;
mod hex;
mod string;

use std::io;

use noodles_sam::alignment::record::data::field::Value;

use self::{array::write_array, hex::write_hex, string::write_string};
use crate::record::codec::encoder::num::{
    write_f32_le, write_i8, write_i16_le, write_i32_le, write_u8, write_u16_le, write_u32_le,
};

pub fn write_value(dst: &mut Vec<u8>, value: &Value) -> io::Result<()> {
    match value {
        Value::Character(c) => write_u8(dst, *c),
        Value::Int8(n) => write_i8(dst, *n),
        Value::UInt8(n) => write_u8(dst, *n),
        Value::Int16(n) => write_i16_le(dst, *n),
        Value::UInt16(n) => write_u16_le(dst, *n),
        Value::Int32(n) => write_i32_le(dst, *n),
        Value::UInt32(n) => write_u32_le(dst, *n),
        Value::Float(n) => write_f32_le(dst, *n),
        Value::String(s) => write_string(dst, s)?,
        Value::Hex(s) => write_hex(dst, s)?,
        Value::Array(array) => write_array(dst, array)?,
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;

    use super::*;

    #[test]
    fn test_write_value() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_sam::alignment::record::data::field::value::Array;

        fn t(buf: &mut Vec<u8>, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Value::Character(b'n'), b"n")?;
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
