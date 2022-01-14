mod subtype;
mod ty;

pub use self::{subtype::read_subtype, ty::read_type};

use std::io::{self, BufRead};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::record::data::field::{
    value::{Subtype, Type},
    Value,
};

/// Reads a BAM record data field value.
///
/// The stream is expected to be at the start of the value, i.e., after the tag and data type.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use noodles_bam::{
///     reader::record::data::field::read_value,
///     record::data::field::{value::Type, Value}
/// };
///
/// let data = [0x01, 0x00, 0x00, 0x00];
/// let mut reader = &data[..];
///
/// assert_eq!(
///     read_value(&mut reader, Type::Int32)?,
///     Value::Int32(1)
/// );
/// # Ok::<(), io::Error>(())
/// ```
pub fn read_value<R>(reader: &mut R, ty: Type) -> io::Result<Value>
where
    R: BufRead,
{
    match ty {
        Type::Char => reader.read_u8().map(char::from).map(Value::Char),
        Type::Int8 => reader.read_i8().map(Value::Int8),
        Type::UInt8 => reader.read_u8().map(Value::UInt8),
        Type::Int16 => reader.read_i16::<LittleEndian>().map(Value::Int16),
        Type::UInt16 => reader.read_u16::<LittleEndian>().map(Value::UInt16),
        Type::Int32 => reader.read_i32::<LittleEndian>().map(Value::Int32),
        Type::UInt32 => reader.read_u32::<LittleEndian>().map(Value::UInt32),
        Type::Float => reader.read_f32::<LittleEndian>().map(Value::Float),
        Type::String => read_string(reader).map(Value::String),
        Type::Hex => read_string(reader).map(Value::Hex),
        Type::Array => read_array(reader),
    }
}

fn read_string<R>(reader: &mut R) -> io::Result<String>
where
    R: BufRead,
{
    let mut buf = Vec::new();
    reader.read_until(b'\0', &mut buf)?;
    buf.pop();
    String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_array<R>(reader: &mut R) -> io::Result<Value>
where
    R: BufRead,
{
    let subtype = read_subtype(reader)?;

    let len = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    match subtype {
        Subtype::Int8 => {
            let mut buf = vec![0; len];
            reader.read_i8_into(&mut buf)?;
            Ok(Value::Int8Array(buf))
        }
        Subtype::UInt8 => {
            let mut buf = vec![0; len];
            reader.read_exact(&mut buf)?;
            Ok(Value::UInt8Array(buf))
        }
        Subtype::Int16 => {
            let mut buf = vec![0; len];
            reader.read_i16_into::<LittleEndian>(&mut buf)?;
            Ok(Value::Int16Array(buf))
        }
        Subtype::UInt16 => {
            let mut buf = vec![0; len];
            reader.read_u16_into::<LittleEndian>(&mut buf)?;
            Ok(Value::UInt16Array(buf))
        }
        Subtype::Int32 => {
            let mut buf = vec![0; len];
            reader.read_i32_into::<LittleEndian>(&mut buf)?;
            Ok(Value::Int32Array(buf))
        }
        Subtype::UInt32 => {
            let mut buf = vec![0; len];
            reader.read_u32_into::<LittleEndian>(&mut buf)?;
            Ok(Value::UInt32Array(buf))
        }
        Subtype::Float => {
            let mut buf = vec![0.0; len];
            reader.read_f32_into::<LittleEndian>(&mut buf)?;
            Ok(Value::FloatArray(buf))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_value() -> io::Result<()> {
        fn t(mut data: &[u8], ty: Type, expected: Value) -> io::Result<()> {
            let actual = read_value(&mut data, ty)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[b'n'], Type::Char, Value::Char('n'))?;
        t(&[0x00], Type::Int8, Value::Int8(0))?;
        t(&[0x00], Type::UInt8, Value::UInt8(0))?;
        t(&[0x00, 0x00], Type::Int16, Value::Int16(0))?;
        t(&[0x00, 0x00], Type::UInt16, Value::UInt16(0))?;
        t(&[0x00, 0x00, 0x00, 0x00], Type::Int32, Value::Int32(0))?;
        t(&[0x00, 0x00, 0x00, 0x00], Type::UInt32, Value::UInt32(0))?;
        t(&[0x00, 0x00, 0x00, 0x00], Type::Float, Value::Float(0.0))?;
        t(
            &[0x6e, 0x64, 0x6c, 0x73, 0x00],
            Type::String,
            Value::String(String::from("ndls")),
        )?;
        t(
            &[0x43, 0x41, 0x46, 0x45, 0x00],
            Type::Hex,
            Value::Hex(String::from("CAFE")),
        )?;

        t(
            &[b'c', 0x01, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Int8Array(vec![0]),
        )?;
        t(
            &[b'C', 0x01, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::UInt8Array(vec![0]),
        )?;
        t(
            &[b's', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Int16Array(vec![0]),
        )?;
        t(
            &[b'S', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::UInt16Array(vec![0]),
        )?;
        t(
            &[b'i', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Int32Array(vec![0]),
        )?;
        t(
            &[b'I', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::UInt32Array(vec![0]),
        )?;
        t(
            &[b'f', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::FloatArray(vec![0.0]),
        )?;

        Ok(())
    }
}
