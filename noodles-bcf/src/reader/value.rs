mod ty;

use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use self::ty::{read_type, Type};

#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    Int8(i8),
    Int8Array(Vec<i8>),
    Int16(i16),
    Int16Array(Vec<i16>),
    Int32(i32),
    Int32Array(Vec<i32>),
    Float(f32),
    FloatArray(Vec<f32>),
    String(Option<String>),
}

pub fn read_value<R>(reader: &mut R) -> io::Result<Value>
where
    R: Read,
{
    let ty = read_type(reader)?;

    match ty {
        Type::Int8(len) => match len {
            0 => todo!("unhandled int8 length: {}", len),
            1 => read_i8(reader),
            _ => read_i8_array(reader, len),
        },
        Type::Int16(len) => match len {
            0 => todo!("unhandled int16 length: {}", len),
            1 => read_i16(reader),
            _ => read_i16_array(reader, len),
        },
        Type::Int32(len) => match len {
            0 => todo!("unhandled int32 length: {}", len),
            1 => read_i32(reader),
            _ => read_i32_array(reader, len),
        },
        Type::Float(len) => match len {
            0 => todo!("unhandled float length: {}", len),
            1 => read_float(reader),
            _ => read_float_array(reader, len),
        },
        Type::String(len) => match len {
            0 => Ok(Value::String(None)),
            _ => read_string(reader, len).map(Some).map(Value::String),
        },
    }
}

fn read_i8<R>(reader: &mut R) -> io::Result<Value>
where
    R: Read,
{
    reader.read_i8().map(Value::Int8)
}

fn read_i8_array<R>(reader: &mut R, len: usize) -> io::Result<Value>
where
    R: Read,
{
    let mut buf = vec![0; len];
    reader.read_i8_into(&mut buf)?;
    Ok(Value::Int8Array(buf))
}

fn read_i16<R>(reader: &mut R) -> io::Result<Value>
where
    R: Read,
{
    reader.read_i16::<LittleEndian>().map(Value::Int16)
}

fn read_i16_array<R>(reader: &mut R, len: usize) -> io::Result<Value>
where
    R: Read,
{
    let mut buf = vec![0; len];
    reader.read_i16_into::<LittleEndian>(&mut buf)?;
    Ok(Value::Int16Array(buf))
}

fn read_i32<R>(reader: &mut R) -> io::Result<Value>
where
    R: Read,
{
    reader.read_i32::<LittleEndian>().map(Value::Int32)
}

fn read_i32_array<R>(reader: &mut R, len: usize) -> io::Result<Value>
where
    R: Read,
{
    let mut buf = vec![0; len];
    reader.read_i32_into::<LittleEndian>(&mut buf)?;
    Ok(Value::Int32Array(buf))
}

fn read_float<R>(reader: &mut R) -> io::Result<Value>
where
    R: Read,
{
    reader.read_f32::<LittleEndian>().map(Value::Float)
}

fn read_float_array<R>(reader: &mut R, len: usize) -> io::Result<Value>
where
    R: Read,
{
    let mut buf = vec![0.0; len];
    reader.read_f32_into::<LittleEndian>(&mut buf)?;
    Ok(Value::FloatArray(buf))
}

fn read_string<R>(reader: &mut R, len: usize) -> io::Result<String>
where
    R: Read,
{
    let mut buf = vec![0; len];
    reader.read_exact(&mut buf)?;
    String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_value() {
        let data = [0x11, 0x05];
        let mut reader = &data[..];
        assert!(matches!(read_value(&mut reader), Ok(Value::Int8(5))));

        let data = [0x31, 0x05, 0x08, 0x0d];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Value::Int8Array(values)) if values == vec![5, 8, 13]
        ));

        let data = [0x12, 0x79, 0x01];
        let mut reader = &data[..];
        assert!(matches!(read_value(&mut reader), Ok(Value::Int16(377))));

        let data = [0x32, 0x79, 0x01, 0x62, 0x02, 0xdb, 0x03];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Value::Int16Array(values)) if values == vec![377, 610, 987]
        ));

        let data = [0x13, 0x11, 0x25, 0x01, 0x00];
        let mut reader = &data[..];
        assert!(matches!(read_value(&mut reader), Ok(Value::Int32(75025))));

        let data = [
            0x33, 0x11, 0x25, 0x01, 0x00, 0x31, 0xda, 0x01, 0x00, 0x42, 0xff, 0x02, 0x00,
        ];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Value::Int32Array(values)) if values == vec![75025, 121393, 196418]
        ));

        let data = [0x15, 0x00, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert!(matches!(read_value(&mut reader), Ok(Value::Float(value)) if value == 0.0));

        let data = [0x25, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3f];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Value::FloatArray(value)) if value == vec![0.0, 0.5]
        ));

        let data = [0x07];
        let mut reader = &data[..];
        assert!(matches!(read_value(&mut reader), Ok(Value::String(None))));

        let data = [0x17, 0x6e];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Value::String(Some(value))) if value == "n"
        ));

        let data = [0x47, 0x6e, 0x64, 0x6c, 0x73];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Value::String(Some(value))) if value == "ndls"
        ));
    }
}
