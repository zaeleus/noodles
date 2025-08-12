mod ty;

pub use self::ty::write_type;

use std::io::{self, Write};

use crate::{
    io::writer::num::{write_f32_le, write_i8, write_i16_le, write_i32_le},
    record::codec::value::{Array, Float, Int8, Int16, Int32, Type, Value},
};

pub fn write_value<W>(writer: &mut W, value: Option<Value<'_>>) -> io::Result<()>
where
    W: Write,
{
    match value {
        None => write_missing(writer),
        Some(Value::Int8(v)) => write_int8(writer, v),
        Some(Value::Int16(v)) => write_int16(writer, v),
        Some(Value::Int32(v)) => write_int32(writer, v),
        Some(Value::Float(v)) => write_float(writer, v),
        Some(Value::String(v)) => write_string(writer, v),
        Some(Value::Array(array)) => write_array(writer, &array),
    }
}

fn write_missing<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    write_type(writer, None)
}

fn write_int8<W>(writer: &mut W, value: Option<Int8>) -> io::Result<()>
where
    W: Write,
{
    match value {
        None => {
            write_type(writer, Some(Type::Int8(0)))?;
        }
        Some(v) => {
            write_type(writer, Some(Type::Int8(1)))?;
            write_i8(writer, i8::from(v))?;
        }
    }

    Ok(())
}

fn write_int16<W>(writer: &mut W, value: Option<Int16>) -> io::Result<()>
where
    W: Write,
{
    match value {
        None => {
            write_type(writer, Some(Type::Int16(0)))?;
        }
        Some(v) => {
            write_type(writer, Some(Type::Int16(1)))?;
            write_i16_le(writer, i16::from(v))?;
        }
    }

    Ok(())
}

fn write_int32<W>(writer: &mut W, value: Option<Int32>) -> io::Result<()>
where
    W: Write,
{
    match value {
        None => {
            write_type(writer, Some(Type::Int32(0)))?;
        }
        Some(v) => {
            write_type(writer, Some(Type::Int32(1)))?;
            write_i32_le(writer, i32::from(v))?;
        }
    }

    Ok(())
}

fn write_float<W>(writer: &mut W, value: Option<Float>) -> io::Result<()>
where
    W: Write,
{
    match value {
        None => {
            write_type(writer, Some(Type::Float(0)))?;
        }
        Some(v) => {
            write_type(writer, Some(Type::Float(1)))?;
            write_f32_le(writer, f32::from(v))?;
        }
    }

    Ok(())
}

fn write_string<W>(writer: &mut W, value: Option<&str>) -> io::Result<()>
where
    W: Write,
{
    match value {
        None => {
            write_type(writer, Some(Type::String(0)))?;
        }
        Some(s) => {
            write_type(writer, Some(Type::String(s.len())))?;
            writer.write_all(s.as_bytes())?;
        }
    }

    Ok(())
}

fn write_array<W>(writer: &mut W, array: &Array) -> io::Result<()>
where
    W: Write,
{
    match array {
        Array::Int8(values) => {
            write_type(writer, Some(Type::Int8(values.len())))?;

            for result in values.iter() {
                let n = result?;
                write_i8(writer, n)?;
            }
        }
        Array::Int16(values) => {
            write_type(writer, Some(Type::Int16(values.len())))?;

            for result in values.iter() {
                let n = result?;
                write_i16_le(writer, n)?;
            }
        }
        Array::Int32(values) => {
            write_type(writer, Some(Type::Int32(values.len())))?;

            for result in values.iter() {
                let n = result?;
                write_i32_le(writer, n)?;
            }
        }
        Array::Float(values) => {
            write_type(writer, Some(Type::Float(values.len())))?;

            for result in values.iter() {
                let n = result?;
                write_f32_le(writer, n)?;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_value() -> io::Result<()> {
        let mut buf = Vec::new();
        write_value(&mut buf, None)?;
        assert_eq!(buf, [0x00]);

        buf.clear();
        write_value(&mut buf, Some(Value::Int8(None)))?;
        assert_eq!(buf, [0x01]);

        buf.clear();
        write_value(&mut buf, Some(Value::Int8(Some(Int8::Value(5)))))?;
        assert_eq!(buf, [0x11, 0x05]);

        buf.clear();
        write_value(&mut buf, Some(Value::Int16(None)))?;
        assert_eq!(buf, [0x02]);

        buf.clear();
        write_value(&mut buf, Some(Value::Int16(Some(Int16::Value(377)))))?;
        assert_eq!(buf, [0x12, 0x79, 0x01]);

        buf.clear();
        write_value(&mut buf, Some(Value::Int32(None)))?;
        assert_eq!(buf, [0x03]);

        buf.clear();
        write_value(&mut buf, Some(Value::Int32(Some(Int32::Value(75025)))))?;
        assert_eq!(buf, [0x13, 0x11, 0x25, 0x01, 0x00]);

        buf.clear();
        write_value(&mut buf, Some(Value::Float(Some(Float::Value(0.0)))))?;
        assert_eq!(buf, [0x15, 0x00, 0x00, 0x00, 0x00]);

        buf.clear();
        write_value(&mut buf, Some(Value::String(None)))?;
        assert_eq!(buf, [0x07]);

        buf.clear();
        write_value(&mut buf, Some(Value::String(Some("n"))))?;
        assert_eq!(buf, [0x17, 0x6e]);

        buf.clear();
        write_value(&mut buf, Some(Value::String(Some("ndls"))))?;
        assert_eq!(buf, [0x47, 0x6e, 0x64, 0x6c, 0x73]);

        buf.clear();
        write_value(
            &mut buf,
            Some(Value::Array(Array::Int8(Box::new(vec![5, 8, 13])))),
        )?;
        assert_eq!(buf, [0x31, 0x05, 0x08, 0x0d]);

        buf.clear();
        write_value(
            &mut buf,
            Some(Value::Array(Array::Int16(Box::new(vec![377, 610, 987])))),
        )?;
        assert_eq!(buf, [0x32, 0x79, 0x01, 0x62, 0x02, 0xdb, 0x03]);

        buf.clear();
        write_value(
            &mut buf,
            Some(Value::Array(Array::Int32(Box::new(vec![
                75025, 121393, 196418,
            ])))),
        )?;
        assert_eq!(
            buf,
            [
                0x33, 0x11, 0x25, 0x01, 0x00, 0x31, 0xda, 0x01, 0x00, 0x42, 0xff, 0x02, 0x00
            ]
        );

        buf.clear();
        write_value(
            &mut buf,
            Some(Value::Array(Array::Float(Box::new(vec![0.0, 0.5])))),
        )?;
        assert_eq!(buf, [0x25, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3f]);

        Ok(())
    }
}
