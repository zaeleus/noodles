mod subtype;
mod ty;

pub use self::ty::write_type;

use std::io::{self, Write};

use self::subtype::write_subtype;
use crate::{
    record::data::field::{value::Character, value::Subtype, Value},
    writer::num,
};

pub fn write_value<W>(writer: &mut W, value: &Value) -> io::Result<()>
where
    W: Write,
{
    const ARRAY_VALUE_DELIMITER: u8 = b',';

    match value {
        Value::Character(c) => write_character(writer, *c),
        Value::Int8(n) => num::write_i8(writer, *n),
        Value::UInt8(n) => num::write_u8(writer, *n),
        Value::Int16(n) => num::write_i16(writer, *n),
        Value::UInt16(n) => num::write_u16(writer, *n),
        Value::Int32(n) => num::write_i32(writer, *n),
        Value::UInt32(n) => num::write_u32(writer, *n),
        Value::Float(n) => num::write_f32(writer, *n),
        Value::String(s) | Value::Hex(s) => writer.write_all(s.as_bytes()),
        Value::Int8Array(values) => {
            write_subtype(writer, Subtype::Int8)?;

            for &n in values {
                writer.write_all(&[ARRAY_VALUE_DELIMITER])?;
                num::write_i8(writer, n)?;
            }

            Ok(())
        }
        Value::UInt8Array(values) => {
            write_subtype(writer, Subtype::UInt8)?;

            for &n in values {
                writer.write_all(&[ARRAY_VALUE_DELIMITER])?;
                num::write_u8(writer, n)?;
            }

            Ok(())
        }
        Value::Int16Array(values) => {
            write_subtype(writer, Subtype::Int16)?;

            for &n in values {
                writer.write_all(&[ARRAY_VALUE_DELIMITER])?;
                num::write_i16(writer, n)?;
            }

            Ok(())
        }
        Value::UInt16Array(values) => {
            write_subtype(writer, Subtype::UInt16)?;

            for &n in values {
                writer.write_all(&[ARRAY_VALUE_DELIMITER])?;
                num::write_u16(writer, n)?;
            }

            Ok(())
        }
        Value::Int32Array(values) => {
            write_subtype(writer, Subtype::Int32)?;

            for &n in values {
                writer.write_all(&[ARRAY_VALUE_DELIMITER])?;
                num::write_i32(writer, n)?;
            }

            Ok(())
        }
        Value::UInt32Array(values) => {
            write_subtype(writer, Subtype::UInt32)?;

            for &n in values {
                writer.write_all(&[ARRAY_VALUE_DELIMITER])?;
                num::write_u32(writer, n)?;
            }

            Ok(())
        }
        Value::FloatArray(values) => {
            write_subtype(writer, Subtype::Float)?;

            for &n in values {
                write!(writer, ",{}", n)?;
            }

            Ok(())
        }
    }
}

pub fn write_character<W>(writer: &mut W, character: Character) -> io::Result<()>
where
    W: Write,
{
    let c = u8::from(character);
    writer.write_all(&[c])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_value() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(
            &mut buf,
            &Value::Character(Character::try_from('n')?),
            &[b'n'],
        )?;
        t(&mut buf, &Value::Int8(1), b"1")?;
        t(&mut buf, &Value::UInt8(2), b"2")?;
        t(&mut buf, &Value::Int16(3), b"3")?;
        t(&mut buf, &Value::UInt16(5), b"5")?;
        t(&mut buf, &Value::Int32(8), b"8")?;
        t(&mut buf, &Value::UInt32(13), b"13")?;
        t(&mut buf, &Value::Float(8.0), b"8")?;
        t(&mut buf, &Value::String(String::from("ndls")), b"ndls")?;
        t(&mut buf, &Value::Hex(String::from("CAFE")), b"CAFE")?;
        t(&mut buf, &Value::Int8Array(vec![1, -2]), b"c,1,-2")?;
        t(&mut buf, &Value::UInt8Array(vec![3, 5]), b"C,3,5")?;
        t(&mut buf, &Value::Int16Array(vec![8, -13]), b"s,8,-13")?;
        t(&mut buf, &Value::UInt16Array(vec![21, 34]), b"S,21,34")?;
        t(&mut buf, &Value::Int32Array(vec![55, -89]), b"i,55,-89")?;
        t(&mut buf, &Value::UInt32Array(vec![144, 223]), b"I,144,223")?;
        t(&mut buf, &Value::FloatArray(vec![8.0, 13.0]), b"f,8,13")?;

        Ok(())
    }
}
