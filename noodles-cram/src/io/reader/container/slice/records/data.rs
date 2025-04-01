use std::{io, ops::Range};

use bstr::{BStr, ByteSlice};
use noodles_sam::alignment::record::data::field::Type;

use crate::record::data::field::{
    value::{array::Values, Array},
    Value,
};

pub(super) fn read_value(src: &[u8], ty: Type) -> io::Result<Value<'_>> {
    match ty {
        Type::Character => read_u8(src).map(Value::Character),
        Type::Int8 => read_u8(src).map(|n| Value::Int8(n as i8)),
        Type::UInt8 => read_u8(src).map(Value::UInt8),
        Type::Int16 => read_u16_le(src).map(|n| Value::Int16(n as i16)),
        Type::UInt16 => read_u16_le(src).map(Value::UInt16),
        Type::Int32 => read_u32_le(src).map(|n| Value::Int32(n as i32)),
        Type::UInt32 => read_u32_le(src).map(Value::UInt32),
        Type::Float => read_f32_le(src).map(Value::Float),
        Type::String => read_string(src).map(Value::String),
        Type::Hex => read_string(src).map(Value::Hex),
        Type::Array => read_array(src).map(Value::Array),
    }
}

fn read_u8(src: &[u8]) -> io::Result<u8> {
    src.split_first()
        .map(|(n, _)| n)
        .copied()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
}

fn read_u16_le(src: &[u8]) -> io::Result<u16> {
    src.split_first_chunk()
        .map(|(buf, _)| u16::from_le_bytes(*buf))
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
}

fn read_u32_le(src: &[u8]) -> io::Result<u32> {
    src.split_first_chunk()
        .map(|(buf, _)| u32::from_le_bytes(*buf))
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
}

fn read_f32_le(src: &[u8]) -> io::Result<f32> {
    src.split_first_chunk()
        .map(|(buf, _)| f32::from_le_bytes(*buf))
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
}

fn read_string(src: &[u8]) -> io::Result<&BStr> {
    const NUL: u8 = 0x00;

    src.strip_suffix(&[NUL])
        .map(|s| s.as_bstr())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing NUL terminator"))
}

fn read_array(src: &[u8]) -> io::Result<Array<'_>> {
    const LENGTH_RANGE: Range<usize> = 1..5;

    let subtype = src
        .first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    let buf = src
        .get(LENGTH_RANGE)
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;
    let n = u32::from_le_bytes(buf.try_into().unwrap());
    let len = usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    match subtype {
        b'c' => Ok(Array::Int8(Values::new(src, len))),
        b'C' => Ok(Array::UInt8(Values::new(src, len))),
        b's' => Ok(Array::Int16(Values::new(src, len))),
        b'S' => Ok(Array::UInt16(Values::new(src, len))),
        b'i' => Ok(Array::Int32(Values::new(src, len))),
        b'I' => Ok(Array::UInt32(Values::new(src, len))),
        b'f' => Ok(Array::Float(Values::new(src, len))),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid subtype",
        )),
    }
}
