mod array;
mod character;

use std::io::{self, Write};

use self::{array::write_array, character::write_character};
use crate::{alignment::record::data::field::Value, io::writer::num};

pub fn write_value<W>(writer: &mut W, value: &Value) -> io::Result<()>
where
    W: Write,
{
    match value {
        Value::Character(c) => write_character(writer, *c),
        Value::Int8(n) => num::write_i8(writer, *n),
        Value::UInt8(n) => num::write_u8(writer, *n),
        Value::Int16(n) => num::write_i16(writer, *n),
        Value::UInt16(n) => num::write_u16(writer, *n),
        Value::Int32(n) => num::write_i32(writer, *n),
        Value::UInt32(n) => num::write_u32(writer, *n),
        Value::Float(n) => num::write_f32(writer, *n),
        Value::String(s) => writer.write_all(s),
        Value::Hex(s) => writer.write_all(s),
        Value::Array(array) => write_array(writer, array),
    }
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;

    use super::*;

    #[test]
    fn test_write_value() -> Result<(), Box<dyn std::error::Error>> {
        use super::array::tests::T;
        use crate::alignment::record::data::field::value::Array;

        fn t(buf: &mut Vec<u8>, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Value::Character(b'n'), &[b'n'])?;
        t(&mut buf, &Value::Int8(1), b"1")?;
        t(&mut buf, &Value::UInt8(2), b"2")?;
        t(&mut buf, &Value::Int16(3), b"3")?;
        t(&mut buf, &Value::UInt16(5), b"5")?;
        t(&mut buf, &Value::Int32(8), b"8")?;
        t(&mut buf, &Value::UInt32(13), b"13")?;
        t(&mut buf, &Value::Float(8.0), b"8")?;
        t(&mut buf, &Value::String(b"ndls".as_bstr()), b"ndls")?;
        t(&mut buf, &Value::Hex(b"CAFE".as_bstr()), b"CAFE")?;
        t(
            &mut buf,
            &Value::Array(Array::Int8(Box::new(T::new(&[0])))),
            b"c,0",
        )?;

        Ok(())
    }
}
