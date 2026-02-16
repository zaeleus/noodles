use std::{io, mem, ops::Range};

use bstr::BString;
use noodles_sam::{
    self as sam,
    alignment::{record::data::field::Type, record_buf::data::field::Value as ValueBuf},
};

pub(super) fn read_value_buf(src: &[u8], ty: Type) -> io::Result<ValueBuf> {
    match ty {
        Type::Character => read_u8(src).map(ValueBuf::Character),
        Type::Int8 => read_u8(src).map(|n| ValueBuf::Int8(n as i8)),
        Type::UInt8 => read_u8(src).map(ValueBuf::UInt8),
        Type::Int16 => read_u16_le(src).map(|n| ValueBuf::Int16(n as i16)),
        Type::UInt16 => read_u16_le(src).map(ValueBuf::UInt16),
        Type::Int32 => read_u32_le(src).map(|n| ValueBuf::Int32(n as i32)),
        Type::UInt32 => read_u32_le(src).map(ValueBuf::UInt32),
        Type::Float => read_f32_le(src).map(ValueBuf::Float),
        Type::String => read_string_owned(src).map(ValueBuf::String),
        Type::Hex => read_string_owned(src).map(ValueBuf::Hex),
        Type::Array => read_array_buf(src),
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

fn read_string_owned(src: &[u8]) -> io::Result<BString> {
    const NUL: u8 = 0x00;

    src.strip_suffix(&[NUL])
        .map(BString::from)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing NUL terminator"))
}

fn read_array_buf(src: &[u8]) -> io::Result<ValueBuf> {
    use sam::alignment::record_buf::data::field::value::Array as ArrayBuf;

    const OFFSET: usize = 5;
    const LENGTH_RANGE: Range<usize> = 1..5;

    let subtype = src
        .first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    let buf = src
        .get(LENGTH_RANGE)
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;
    let n = u32::from_le_bytes(buf.try_into().unwrap());
    let len = usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let data = &src[OFFSET..];

    let element_size = match subtype {
        b'c' | b'C' => 1,
        b's' | b'S' => mem::size_of::<i16>(),
        b'i' | b'I' => mem::size_of::<i32>(),
        b'f' => mem::size_of::<f32>(),
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid subtype",
            ));
        }
    };

    let expected = len
        .checked_mul(element_size)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "array length overflow"))?;

    if data.len() < expected {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    match subtype {
        b'c' => {
            let values: Vec<i8> = data[..len].iter().map(|&b| b as i8).collect();
            Ok(ValueBuf::Array(ArrayBuf::Int8(values)))
        }
        b'C' => {
            let values: Vec<u8> = data[..len].to_vec();
            Ok(ValueBuf::Array(ArrayBuf::UInt8(values)))
        }
        b's' => {
            let values: Vec<i16> = data[..expected]
                .chunks_exact(mem::size_of::<i16>())
                .map(|buf| i16::from_le_bytes(buf.try_into().unwrap()))
                .collect();
            Ok(ValueBuf::Array(ArrayBuf::Int16(values)))
        }
        b'S' => {
            let values: Vec<u16> = data[..expected]
                .chunks_exact(mem::size_of::<u16>())
                .map(|buf| u16::from_le_bytes(buf.try_into().unwrap()))
                .collect();
            Ok(ValueBuf::Array(ArrayBuf::UInt16(values)))
        }
        b'i' => {
            let values: Vec<i32> = data[..expected]
                .chunks_exact(mem::size_of::<i32>())
                .map(|buf| i32::from_le_bytes(buf.try_into().unwrap()))
                .collect();
            Ok(ValueBuf::Array(ArrayBuf::Int32(values)))
        }
        b'I' => {
            let values: Vec<u32> = data[..expected]
                .chunks_exact(mem::size_of::<u32>())
                .map(|buf| u32::from_le_bytes(buf.try_into().unwrap()))
                .collect();
            Ok(ValueBuf::Array(ArrayBuf::UInt32(values)))
        }
        b'f' => {
            let values: Vec<f32> = data[..expected]
                .chunks_exact(mem::size_of::<f32>())
                .map(|buf| f32::from_le_bytes(buf.try_into().unwrap()))
                .collect();
            Ok(ValueBuf::Array(ArrayBuf::Float(values)))
        }
        _ => unreachable!(),
    }
}
