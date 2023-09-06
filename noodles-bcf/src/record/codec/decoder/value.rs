mod ty;

pub use self::ty::read_type;

use std::{io, str};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::lazy::record::{
    value::{Array, Float, Int16, Int32, Int8, Type},
    Value,
};

pub fn read_value<'a>(src: &mut &'a [u8]) -> io::Result<Option<Value<'a>>> {
    let ty = read_type(src)?;

    match ty {
        Some(Type::Int8(len)) => match len {
            0 => Ok(Some(Value::Int8(None))),
            1 => read_i8(src)
                .map(Int8::from)
                .map(Some)
                .map(Value::Int8)
                .map(Some),
            _ => read_i8_array(src, len)
                .map(Array::Int8)
                .map(Value::Array)
                .map(Some),
        },
        Some(Type::Int16(len)) => match len {
            0 => Ok(Some(Value::Int16(None))),
            1 => read_i16(src)
                .map(Int16::from)
                .map(Some)
                .map(Value::Int16)
                .map(Some),
            _ => read_i16_array(src, len)
                .map(Array::Int16)
                .map(Value::Array)
                .map(Some),
        },
        Some(Type::Int32(len)) => match len {
            0 => Ok(Some(Value::Int32(None))),
            1 => read_i32(src)
                .map(Int32::from)
                .map(Some)
                .map(Value::Int32)
                .map(Some),
            _ => read_i32_array(src, len)
                .map(Array::Int32)
                .map(Value::Array)
                .map(Some),
        },
        Some(Type::Float(len)) => match len {
            0 => Ok(Some(Value::Float(None))),
            1 => read_float(src)
                .map(Float::from)
                .map(Some)
                .map(Value::Float)
                .map(Some),
            _ => read_float_array(src, len)
                .map(Array::Float)
                .map(Value::Array)
                .map(Some),
        },
        Some(Type::String(len)) => match len {
            0 => Ok(Some(Value::String(None))),
            _ => read_string(src, len).map(Some).map(Value::String).map(Some),
        },
        None => Ok(None),
    }
}

fn read_i8(src: &mut &[u8]) -> io::Result<i8> {
    src.read_i8()
}

fn read_i8_array(src: &mut &[u8], len: usize) -> io::Result<Vec<i8>> {
    let mut buf = vec![0; len];
    src.read_i8_into(&mut buf)?;
    Ok(buf)
}

fn read_i16(src: &mut &[u8]) -> io::Result<i16> {
    src.read_i16::<LittleEndian>()
}

fn read_i16_array(src: &mut &[u8], len: usize) -> io::Result<Vec<i16>> {
    let mut buf = vec![0; len];
    src.read_i16_into::<LittleEndian>(&mut buf)?;
    Ok(buf)
}

fn read_i32(src: &mut &[u8]) -> io::Result<i32> {
    src.read_i32::<LittleEndian>()
}

fn read_i32_array(src: &mut &[u8], len: usize) -> io::Result<Vec<i32>> {
    let mut buf = vec![0; len];
    src.read_i32_into::<LittleEndian>(&mut buf)?;
    Ok(buf)
}

fn read_float(src: &mut &[u8]) -> io::Result<f32> {
    src.read_f32::<LittleEndian>()
}

fn read_float_array(src: &mut &[u8], len: usize) -> io::Result<Vec<f32>> {
    let mut buf = vec![0.0; len];
    src.read_f32_into::<LittleEndian>(&mut buf)?;
    Ok(buf)
}

fn read_string<'a>(src: &mut &'a [u8], len: usize) -> io::Result<&'a str> {
    if src.len() < len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let (buf, rest) = src.split_at(len);
    let s = str::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    *src = rest;

    Ok(s)
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
            Ok(Some(Value::Array(Array::Int8(values)))) if values == [5, 8, 13]
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
            Ok(Some(Value::Array(Array::Int16(values)))) if values == [377, 610, 987]
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
            Ok(Some(Value::Array(Array::Int32(values)))) if values == [75025, 121393, 196418]
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
            Ok(Some(Value::Array(Array::Float(value)))) if value == [0.0, 0.5]
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
