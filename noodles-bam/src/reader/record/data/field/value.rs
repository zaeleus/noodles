mod subtype;
mod ty;

pub use self::{subtype::get_subtype, ty::get_type};

use std::{io, mem};

use bytes::Buf;
use noodles_sam::record::data::field::{
    value::Character,
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
/// use noodles_bam::reader::record::data::field::get_value;
/// use noodles_sam::record::data::field::{value::Type, Value};
///
/// let data = [0x01, 0x00, 0x00, 0x00];
/// let mut reader = &data[..];
///
/// assert_eq!(get_value(&mut reader, Type::Int32)?, Value::Int32(1));
/// # Ok::<(), io::Error>(())
/// ```
pub fn get_value<B>(src: &mut B, ty: Type) -> io::Result<Value>
where
    B: Buf,
{
    match ty {
        Type::Character => get_char_value(src),
        Type::Int8 => get_i8_value(src),
        Type::UInt8 => get_u8_value(src),
        Type::Int16 => get_i16_value(src),
        Type::UInt16 => get_u16_value(src),
        Type::Int32 => get_i32_value(src),
        Type::UInt32 => get_u32_value(src),
        Type::Float => get_f32_value(src),
        Type::String => get_string(src).map(Value::String),
        Type::Hex => get_string(src).map(Value::Hex),
        Type::Array => get_array_value(src),
    }
}

fn get_char_value<B>(src: &mut B) -> io::Result<Value>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u8>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Character::try_from(src.get_u8())
        .map(Value::Character)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn get_i8_value<B>(src: &mut B) -> io::Result<Value>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<i8>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(Value::Int8(src.get_i8()))
}

fn get_u8_value<B>(src: &mut B) -> io::Result<Value>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u8>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(Value::UInt8(src.get_u8()))
}

fn get_i16_value<B>(src: &mut B) -> io::Result<Value>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<i16>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(Value::Int16(src.get_i16_le()))
}

fn get_u16_value<B>(src: &mut B) -> io::Result<Value>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u16>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(Value::UInt16(src.get_u16_le()))
}

fn get_i32_value<B>(src: &mut B) -> io::Result<Value>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<i32>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(Value::Int32(src.get_i32_le()))
}

fn get_u32_value<B>(src: &mut B) -> io::Result<Value>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u32>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(Value::UInt32(src.get_u32_le()))
}

fn get_f32_value<B>(src: &mut B) -> io::Result<Value>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<f32>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(Value::Float(src.get_f32_le()))
}

fn get_string<B>(src: &mut B) -> io::Result<String>
where
    B: Buf,
{
    const NUL: u8 = 0x00;

    let len = src.chunk().iter().position(|&b| b == NUL).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "string value missing NUL terminator",
        )
    })?;

    let mut buf = vec![0; len];
    src.copy_to_slice(&mut buf);
    src.advance(1); // Discard the NUL terminator.

    String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn get_array_value<B>(src: &mut B) -> io::Result<Value>
where
    B: Buf,
{
    let subtype = get_subtype(src)?;

    let len = usize::try_from(src.get_i32_le())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    match subtype {
        Subtype::Int8 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_i8());
            }

            Ok(Value::Int8Array(buf))
        }
        Subtype::UInt8 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_u8());
            }

            Ok(Value::UInt8Array(buf))
        }
        Subtype::Int16 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_i16_le());
            }

            Ok(Value::Int16Array(buf))
        }
        Subtype::UInt16 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_u16_le());
            }

            Ok(Value::UInt16Array(buf))
        }
        Subtype::Int32 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_i32_le());
            }

            Ok(Value::Int32Array(buf))
        }
        Subtype::UInt32 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_u32_le());
            }

            Ok(Value::UInt32Array(buf))
        }
        Subtype::Float => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                buf.push(src.get_f32_le());
            }

            Ok(Value::FloatArray(buf))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_value() -> Result<(), Box<dyn std::error::Error>> {
        fn t(mut data: &[u8], ty: Type, expected: Value) -> io::Result<()> {
            let actual = get_value(&mut data, ty)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(
            &[b'n'],
            Type::Character,
            Value::Character(Character::try_from('n')?),
        )?;
        t(&[0x00], Type::Int8, Value::Int8(0))?;
        t(&[0x00], Type::UInt8, Value::UInt8(0))?;
        t(&[0x00, 0x00], Type::Int16, Value::Int16(0))?;
        t(&[0x00, 0x00], Type::UInt16, Value::UInt16(0))?;
        t(&[0x00, 0x00, 0x00, 0x00], Type::Int32, Value::Int32(0))?;
        t(&[0x00, 0x00, 0x00, 0x00], Type::UInt32, Value::UInt32(0))?;
        t(&[0x00, 0x00, 0x00, 0x00], Type::Float, Value::Float(0.0))?;
        t(
            &[b'n', b'd', b'l', b's', 0x00],
            Type::String,
            Value::String(String::from("ndls")),
        )?;
        t(
            &[b'C', b'A', b'F', b'E', 0x00],
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
