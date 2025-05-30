use std::{io, mem, str};

pub(crate) mod array;
mod ty;

use self::array::Values;
pub(crate) use self::{array::Array, ty::Type, ty::read_type};
use super::codec::value::{Float, Int8, Int16, Int32};

pub(crate) enum Value<'a> {
    Int8(Option<Int8>),
    Int16(Option<Int16>),
    Int32(Option<Int32>),
    Float(Option<Float>),
    String(Option<&'a str>),
    Array(Array<'a>),
}

impl Value<'_> {
    pub(crate) fn as_int(&self) -> Option<i32> {
        match self {
            Self::Int8(Some(Int8::Value(n))) => Some(i32::from(*n)),
            Self::Int16(Some(Int16::Value(n))) => Some(i32::from(*n)),
            Self::Int32(Some(Int32::Value(n))) => Some(*n),
            _ => None,
        }
    }
}

pub(crate) fn read_value<'a>(src: &mut &'a [u8]) -> io::Result<Option<Value<'a>>> {
    match read_type(src)? {
        None => Ok(None),
        Some(Type::Int8(0)) => Ok(Some(Value::Int8(None))),
        Some(Type::Int8(1)) => read_int8_value(src),
        Some(Type::Int8(n)) => read_int8_array_value(src, n),
        Some(Type::Int16(0)) => Ok(Some(Value::Int16(None))),
        Some(Type::Int16(1)) => read_int16_value(src),
        Some(Type::Int16(n)) => read_int16_array_value(src, n),
        Some(Type::Int32(0)) => Ok(Some(Value::Int32(None))),
        Some(Type::Int32(1)) => read_int32_value(src),
        Some(Type::Int32(n)) => read_int32_array_value(src, n),
        Some(Type::Float(0)) => Ok(Some(Value::Float(None))),
        Some(Type::Float(1)) => read_float_value(src),
        Some(Type::Float(n)) => read_float_array_value(src, n),
        Some(Type::String(0)) => Ok(Some(Value::String(None))),
        Some(Type::String(n)) => read_string_value(src, n),
    }
}

fn read_int8_value<'a>(src: &mut &'a [u8]) -> io::Result<Option<Value<'a>>> {
    read_i8(src).map(|n| Some(Value::Int8(Some(Int8::from(n)))))
}

fn read_int16_value<'a>(src: &mut &'a [u8]) -> io::Result<Option<Value<'a>>> {
    read_i16(src).map(|n| Some(Value::Int16(Some(Int16::from(n)))))
}

fn read_int32_value<'a>(src: &mut &'a [u8]) -> io::Result<Option<Value<'a>>> {
    read_i32(src).map(|n| Some(Value::Int32(Some(Int32::from(n)))))
}

fn read_float_value<'a>(src: &mut &'a [u8]) -> io::Result<Option<Value<'a>>> {
    read_f32(src).map(|n| Some(Value::Float(Some(Float::from(n)))))
}

fn read_string_value<'a>(src: &mut &'a [u8], len: usize) -> io::Result<Option<Value<'a>>> {
    let buf = split_to(src, len)?;
    let s = str::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    Ok(Some(Value::String(Some(s))))
}

fn read_int8_array_value<'a>(src: &mut &'a [u8], n: usize) -> io::Result<Option<Value<'a>>> {
    let buf = split_to(src, mem::size_of::<i8>() * n)?;
    Ok(Some(Value::Array(Array::Int8(Values::new(buf)))))
}

fn read_int16_array_value<'a>(src: &mut &'a [u8], n: usize) -> io::Result<Option<Value<'a>>> {
    let buf = split_to(src, mem::size_of::<i16>() * n)?;
    Ok(Some(Value::Array(Array::Int16(Values::new(buf)))))
}

fn read_int32_array_value<'a>(src: &mut &'a [u8], n: usize) -> io::Result<Option<Value<'a>>> {
    let buf = split_to(src, mem::size_of::<i32>() * n)?;
    Ok(Some(Value::Array(Array::Int32(Values::new(buf)))))
}

fn read_float_array_value<'a>(src: &mut &'a [u8], n: usize) -> io::Result<Option<Value<'a>>> {
    let buf = split_to(src, mem::size_of::<f32>() * n)?;
    Ok(Some(Value::Array(Array::Float(Values::new(buf)))))
}

fn read_i8(src: &mut &[u8]) -> io::Result<i8> {
    if let Some((b, rest)) = src.split_first() {
        *src = rest;
        Ok(*b as i8)
    } else {
        Err(io::Error::from(io::ErrorKind::UnexpectedEof))
    }
}

fn read_i16(src: &mut &[u8]) -> io::Result<i16> {
    let buf = split_to(src, mem::size_of::<i16>())?;
    // SAFETY: `buf` is 2 bytes.
    Ok(i16::from_le_bytes(buf.try_into().unwrap()))
}

fn read_i32(src: &mut &[u8]) -> io::Result<i32> {
    let buf = split_to(src, mem::size_of::<i32>())?;
    // SAFETY: `buf` is 4 bytes.
    Ok(i32::from_le_bytes(buf.try_into().unwrap()))
}

fn read_f32(src: &mut &[u8]) -> io::Result<f32> {
    let buf = split_to(src, mem::size_of::<f32>())?;
    // SAFETY: `buf` is 4 bytes.
    Ok(f32::from_le_bytes(buf.try_into().unwrap()))
}

fn split_to<'a>(src: &mut &'a [u8], i: usize) -> io::Result<&'a [u8]> {
    if src.len() < i {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let (buf, rest) = src.split_at(i);
    *src = rest;

    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_value() {
        let mut src = &[0x00][..];
        assert!(matches!(read_value(&mut src), Ok(None)));

        let mut src = &[0x01][..];
        assert!(matches!(read_value(&mut src), Ok(Some(Value::Int8(None)))));

        let mut src = &[0x11, 0x05][..];
        assert!(matches!(
            read_value(&mut src),
            Ok(Some(Value::Int8(Some(Int8::Value(5)))))
        ));

        let mut src = &[0x31, 0x05, 0x08, 0x0d][..];
        match read_value(&mut src) {
            Ok(Some(Value::Array(Array::Int8(values)))) => {
                assert_eq!(
                    values.iter().collect::<Vec<_>>(),
                    [Int8::from(5), Int8::from(8), Int8::from(13)]
                );
            }
            _ => panic!(),
        }

        let mut src = &[0x02][..];
        assert!(matches!(read_value(&mut src), Ok(Some(Value::Int16(None)))));

        let mut src = &[0x12, 0x79, 0x01][..];
        assert!(matches!(
            read_value(&mut src),
            Ok(Some(Value::Int16(Some(Int16::Value(377)))))
        ));

        let mut src = &[0x32, 0x79, 0x01, 0x62, 0x02, 0xdb, 0x03][..];
        match read_value(&mut src) {
            Ok(Some(Value::Array(Array::Int16(values)))) => {
                assert_eq!(
                    values.iter().collect::<Vec<_>>(),
                    [Int16::from(377), Int16::from(610), Int16::from(987)]
                );
            }
            _ => panic!(),
        }

        let mut src = &[0x03][..];
        assert!(matches!(read_value(&mut src), Ok(Some(Value::Int32(None)))));

        let mut src = &[0x13, 0x11, 0x25, 0x01, 0x00][..];
        assert!(matches!(
            read_value(&mut src),
            Ok(Some(Value::Int32(Some(Int32::Value(75025)))))
        ));

        let mut src = &[
            0x33, 0x11, 0x25, 0x01, 0x00, 0x31, 0xda, 0x01, 0x00, 0x42, 0xff, 0x02, 0x00,
        ][..];
        match read_value(&mut src) {
            Ok(Some(Value::Array(Array::Int32(values)))) => {
                assert_eq!(
                    values.iter().collect::<Vec<_>>(),
                    [Int32::from(75025), Int32::from(121393), Int32::from(196418)]
                );
            }
            _ => panic!(),
        }

        let mut src = &[0x05][..];
        assert!(matches!(read_value(&mut src), Ok(Some(Value::Float(None)))));

        let mut src = &[0x15, 0x00, 0x00, 0x00, 0x00][..];
        assert!(matches!(
            read_value(&mut src),
            Ok(Some(Value::Float(Some(n)))) if n == Float::from(0.0)
        ));

        let mut src = &[0x25, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3f][..];
        match read_value(&mut src) {
            Ok(Some(Value::Array(Array::Float(values)))) => {
                assert_eq!(
                    values.iter().collect::<Vec<_>>(),
                    [Float::from(0.0), Float::from(0.5)]
                );
            }
            _ => panic!(),
        }

        let mut src = &[0x07][..];
        assert!(matches!(
            read_value(&mut src),
            Ok(Some(Value::String(None)))
        ));

        let mut src = &[0x17, b'n'][..];
        assert!(matches!(
            read_value(&mut src),
            Ok(Some(Value::String(Some("n"))))
        ));

        let mut src = &[0x47, b'n', b'd', b'l', b's'][..];
        assert!(matches!(
            read_value(&mut src),
            Ok(Some(Value::String(Some("ndls"))))
        ));
    }
}
