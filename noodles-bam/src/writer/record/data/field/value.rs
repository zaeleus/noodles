mod subtype;
mod ty;

pub use self::{subtype::write_subtype, ty::write_type};

use std::{
    ffi::CString,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_sam::record::data::field::{value::Subtype, Value};

pub fn write_value<W>(writer: &mut W, value: &Value) -> io::Result<()>
where
    W: Write,
{
    match value {
        Value::Char(c) => writer.write_u8(*c as u8)?,
        Value::Int8(n) => writer.write_i8(*n)?,
        Value::UInt8(n) => writer.write_u8(*n)?,
        Value::Int16(n) => writer.write_i16::<LittleEndian>(*n)?,
        Value::UInt16(n) => writer.write_u16::<LittleEndian>(*n)?,
        Value::Int32(n) => writer.write_i32::<LittleEndian>(*n)?,
        Value::UInt32(n) => writer.write_u32::<LittleEndian>(*n)?,
        Value::Float(n) => writer.write_f32::<LittleEndian>(*n)?,
        Value::String(s) | Value::Hex(s) => write_string(writer, s)?,
        Value::Int8Array(values) => {
            write_array_header(writer, Subtype::Int8, values.len())?;

            for &n in values {
                writer.write_i8(n)?;
            }
        }
        Value::UInt8Array(values) => {
            write_array_header(writer, Subtype::UInt8, values.len())?;

            for &n in values {
                writer.write_u8(n)?;
            }
        }
        Value::Int16Array(values) => {
            write_array_header(writer, Subtype::Int16, values.len())?;

            for &n in values {
                writer.write_i16::<LittleEndian>(n)?;
            }
        }
        Value::UInt16Array(values) => {
            write_array_header(writer, Subtype::UInt16, values.len())?;

            for &n in values {
                writer.write_u16::<LittleEndian>(n)?;
            }
        }
        Value::Int32Array(values) => {
            write_array_header(writer, Subtype::Int32, values.len())?;

            for &n in values {
                writer.write_i32::<LittleEndian>(n)?;
            }
        }
        Value::UInt32Array(values) => {
            write_array_header(writer, Subtype::UInt32, values.len())?;

            for &n in values {
                writer.write_u32::<LittleEndian>(n)?;
            }
        }
        Value::FloatArray(values) => {
            write_array_header(writer, Subtype::Float, values.len())?;

            for &n in values {
                writer.write_f32::<LittleEndian>(n)?;
            }
        }
    }

    Ok(())
}

fn write_string<W>(writer: &mut W, s: &str) -> io::Result<()>
where
    W: Write,
{
    let c_str = CString::new(s).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_all(c_str.as_bytes_with_nul())
}

fn write_array_header<W>(writer: &mut W, subtype: Subtype, len: usize) -> io::Result<()>
where
    W: Write,
{
    write_subtype(writer, subtype)?;

    let n = u32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(n)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_value() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Value::Char('n'), &[b'n'])?;
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
            &Value::Hex(String::from("CAFE")),
            &[b'C', b'A', b'F', b'E', 0x00],
        )?;

        t(
            &mut buf,
            &Value::Int8Array(vec![1, -2]),
            &[
                b'c', // subtype = Int8
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x01, // values[0] = 1
                0xfe, // values[1] = -2
            ],
        )?;

        t(
            &mut buf,
            &Value::UInt8Array(vec![3, 5]),
            &[
                b'C', // subtype = UInt8
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x03, // values[0] = 3
                0x05, // values[1] = 5
            ],
        )?;

        t(
            &mut buf,
            &Value::Int16Array(vec![8, -13]),
            &[
                b's', // subtype = Int16
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x08, 0x00, // values[0] = 8
                0xf3, 0xff, // values[1] = -13
            ],
        )?;

        t(
            &mut buf,
            &Value::UInt16Array(vec![21, 34]),
            &[
                b'S', // subtype = UInt16
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x15, 0x00, // values[0] = 21
                0x22, 0x00, // values[1] = 34
            ],
        )?;

        t(
            &mut buf,
            &Value::Int32Array(vec![55, -89]),
            &[
                b'i', // subtype = Int32
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x37, 0x00, 0x00, 0x00, // values[0] = 55
                0xa7, 0xff, 0xff, 0xff, // values[1] = -89
            ],
        )?;

        t(
            &mut buf,
            &Value::UInt32Array(vec![144, 223]),
            &[
                b'I', // subtype = UInt32
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x90, 0x00, 0x00, 0x00, // values[0] = 55
                0xdf, 0x00, 0x00, 0x00, // values[1] = -89
            ],
        )?;

        t(
            &mut buf,
            &Value::FloatArray(vec![8.0, 13.0]),
            &[
                b'f', // subtype = Float
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x00, 0x00, 0x00, 0x41, // values[0] = 8.0
                0x00, 0x00, 0x50, 0x41, // values[1] = 13.0
            ],
        )?;

        Ok(())
    }
}
