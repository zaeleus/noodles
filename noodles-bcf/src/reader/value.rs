mod ty;

pub use self::ty::read_type;

use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::lazy::record::{
    value::{Float, Int16, Int32, Int8, Type},
    Value,
};

pub fn read_value<R>(reader: &mut R) -> io::Result<Option<Value>>
where
    R: Read,
{
    let ty = read_type(reader)?;

    match ty {
        Some(Type::Int8(len)) => match len {
            0 => Ok(Some(Value::Int8(None))),
            1 => read_i8(reader)
                .map(Int8::from)
                .map(Some)
                .map(Value::Int8)
                .map(Some),
            _ => read_i8_array(reader, len).map(Value::Int8Array).map(Some),
        },
        Some(Type::Int16(len)) => match len {
            0 => Ok(Some(Value::Int16(None))),
            1 => read_i16(reader)
                .map(Int16::from)
                .map(Some)
                .map(Value::Int16)
                .map(Some),
            _ => read_i16_array(reader, len).map(Value::Int16Array).map(Some),
        },
        Some(Type::Int32(len)) => match len {
            0 => Ok(Some(Value::Int32(None))),
            1 => read_i32(reader)
                .map(Int32::from)
                .map(Some)
                .map(Value::Int32)
                .map(Some),
            _ => read_i32_array(reader, len).map(Value::Int32Array).map(Some),
        },
        Some(Type::Float(len)) => match len {
            0 => Ok(Some(Value::Float(None))),
            1 => read_float(reader)
                .map(Float::from)
                .map(Some)
                .map(Value::Float)
                .map(Some),
            _ => read_float_array(reader, len)
                .map(Value::FloatArray)
                .map(Some),
        },
        Some(Type::String(len)) => match len {
            0 => Ok(Some(Value::String(None))),
            _ => read_string(reader, len)
                .map(Some)
                .map(Value::String)
                .map(Some),
        },
        None => Ok(None),
    }
}

fn read_i8<R>(reader: &mut R) -> io::Result<i8>
where
    R: Read,
{
    reader.read_i8()
}

fn read_i8_array<R>(reader: &mut R, len: usize) -> io::Result<Vec<i8>>
where
    R: Read,
{
    let mut buf = vec![0; len];
    reader.read_i8_into(&mut buf)?;
    Ok(buf)
}

fn read_i16<R>(reader: &mut R) -> io::Result<i16>
where
    R: Read,
{
    reader.read_i16::<LittleEndian>()
}

fn read_i16_array<R>(reader: &mut R, len: usize) -> io::Result<Vec<i16>>
where
    R: Read,
{
    let mut buf = vec![0; len];
    reader.read_i16_into::<LittleEndian>(&mut buf)?;
    Ok(buf)
}

fn read_i32<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    reader.read_i32::<LittleEndian>()
}

fn read_i32_array<R>(reader: &mut R, len: usize) -> io::Result<Vec<i32>>
where
    R: Read,
{
    let mut buf = vec![0; len];
    reader.read_i32_into::<LittleEndian>(&mut buf)?;
    Ok(buf)
}

fn read_float<R>(reader: &mut R) -> io::Result<f32>
where
    R: Read,
{
    reader.read_f32::<LittleEndian>()
}

fn read_float_array<R>(reader: &mut R, len: usize) -> io::Result<Vec<f32>>
where
    R: Read,
{
    let mut buf = vec![0.0; len];
    reader.read_f32_into::<LittleEndian>(&mut buf)?;
    Ok(buf)
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
        let data = [0x00];
        let mut reader = &data[..];
        assert!(matches!(read_value(&mut reader), Ok(None)));

        let data = [0x01];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::Int8(None)))
        ));

        let data = [0x11, 0x05];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::Int8(Some(Int8::Value(5)))))
        ));

        let data = [0x31, 0x05, 0x08, 0x0d];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::Int8Array(values))) if values == [5, 8, 13]
        ));

        let data = [0x02];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::Int16(None)))
        ));

        let data = [0x12, 0x79, 0x01];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::Int16(Some(Int16::Value(377)))))
        ));

        let data = [0x32, 0x79, 0x01, 0x62, 0x02, 0xdb, 0x03];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::Int16Array(values))) if values == [377, 610, 987]
        ));

        let data = [0x03];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::Int32(None)))
        ));

        let data = [0x13, 0x11, 0x25, 0x01, 0x00];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::Int32(Some(Int32::Value(75025)))))
        ));

        let data = [
            0x33, 0x11, 0x25, 0x01, 0x00, 0x31, 0xda, 0x01, 0x00, 0x42, 0xff, 0x02, 0x00,
        ];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::Int32Array(values))) if values == [75025, 121393, 196418]
        ));

        let data = [0x05];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::Float(None)))
        ));

        let data = [0x15, 0x00, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::Float(Some(Float::Value(value))))) if value == 0.0
        ));

        let data = [0x25, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3f];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::FloatArray(value))) if value == [0.0, 0.5]
        ));

        let data = [0x07];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::String(None)))
        ));

        let data = [0x17, 0x6e];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::String(Some(value)))) if value == "n"
        ));

        let data = [0x47, 0x6e, 0x64, 0x6c, 0x73];
        let mut reader = &data[..];
        assert!(matches!(
            read_value(&mut reader),
            Ok(Some(Value::String(Some(value)))) if value == "ndls"
        ));
    }
}
