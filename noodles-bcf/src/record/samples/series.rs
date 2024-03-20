//! BCF record samples series.

pub mod value;

use std::{io, mem, ops::Range, str};

use noodles_vcf::{
    self as vcf,
    variant::record::samples::series::{value::Array, Value},
};

use crate::record::value::{array::Values, read_type, read_value, Type};

/// A BCF record samples series.
pub struct Series<'r> {
    id: usize,
    ty: Type,
    src: &'r [u8],
}

impl<'r> Series<'r> {
    /// Returns the name.
    pub fn name<'h>(&self, header: &'h vcf::Header) -> io::Result<&'h str> {
        header
            .string_maps()
            .strings()
            .get_index(self.id)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid string map ID"))
    }

    fn len(&self) -> usize {
        match self.ty {
            Type::Int8(len) => len,
            Type::Int16(len) => len,
            Type::Int32(len) => len,
            Type::Float(len) => len,
            Type::String(len) => len,
        }
    }

    /// Returns the value at the given index.
    pub fn get(&self, i: usize) -> Option<Option<Value<'r>>> {
        match self.ty {
            Type::Int8(len) => get_int8_value(self.src, len, i),
            Type::Int16(len) => get_int16_value(self.src, len, i),
            Type::Int32(len) => get_int32_value(self.src, len, i),
            Type::Float(len) => get_float_value(self.src, len, i),
            Type::String(len) => get_string_value(self.src, len, i),
        }
    }
}

impl<'r> vcf::variant::record::samples::Series for Series<'r> {
    fn name<'a, 'h: 'a>(&'a self, header: &'h vcf::Header) -> io::Result<&'a str> {
        header
            .string_maps()
            .strings()
            .get_index(self.id)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid string map ID"))
    }

    fn get<'a, 'h: 'a>(
        &'a self,
        _: &'h vcf::Header,
        i: usize,
    ) -> Option<Option<io::Result<Value<'a>>>> {
        let value = match self.ty {
            Type::Int8(len) => get_int8_value(self.src, len, i),
            Type::Int16(len) => get_int16_value(self.src, len, i),
            Type::Int32(len) => get_int32_value(self.src, len, i),
            Type::Float(len) => get_float_value(self.src, len, i),
            Type::String(len) => get_string_value(self.src, len, i),
        };

        match value {
            Some(Some(value)) => Some(Some(Ok(value))),
            Some(None) => Some(None),
            None => None,
        }
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        _: &'h vcf::Header,
    ) -> Box<dyn Iterator<Item = io::Result<Option<Value<'a>>>> + 'a> {
        Box::new((0..self.len()).map(|i| {
            self.get(i)
                .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))
        }))
    }
}

pub(super) fn read_series<'a>(src: &mut &'a [u8], sample_count: usize) -> io::Result<Series<'a>> {
    fn size_of(ty: Type) -> usize {
        match ty {
            Type::Int8(n) => mem::size_of::<i8>() * n,
            Type::Int16(n) => mem::size_of::<i16>() * n,
            Type::Int32(n) => mem::size_of::<i32>() * n,
            Type::Float(n) => mem::size_of::<f32>() * n,
            Type::String(n) => mem::size_of::<u8>() * n,
        }
    }

    let id = read_string_map_index(src)?;
    let ty = read_type(src)?.expect("invalid type");

    let len = size_of(ty) * sample_count;
    let (buf, rest) = src.split_at(len);

    *src = rest;

    Ok(Series { id, ty, src: buf })
}

fn read_string_map_index(src: &mut &[u8]) -> io::Result<usize> {
    match read_value(src)?.and_then(|v| v.as_int()) {
        Some(i) => usize::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        None => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid string map index",
        )),
    }
}

fn range<N>(i: usize, len: usize) -> Range<usize> {
    let size = mem::size_of::<N>();
    let start = size * i;
    let end = start + size * len;
    start..end
}

fn get_int8_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    use crate::record::codec::value::Int8;

    let src = src.get(range::<i8>(i, len))?;

    let value = if len == 1 {
        match Int8::from(src[0] as i8) {
            Int8::Value(n) => Some(Value::Integer(i32::from(n))),
            Int8::Missing => None,
            Int8::EndOfVector | Int8::Reserved(_) => todo!(),
        }
    } else {
        let values = Values::<'_, i8>::new(src);
        Some(Value::Array(Array::Integer(Box::new(values))))
    };

    Some(value)
}

fn get_int16_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    use crate::record::codec::value::Int16;

    let src = src.get(range::<i16>(i, len))?;

    let value = if len == 1 {
        // SAFETY: `src` is 2 bytes.
        match Int16::from(i16::from_le_bytes(src.try_into().unwrap())) {
            Int16::Value(n) => Some(Value::Integer(i32::from(n))),
            Int16::Missing => None,
            Int16::EndOfVector | Int16::Reserved(_) => todo!(),
        }
    } else {
        let values = Values::<'_, i16>::new(src);
        Some(Value::Array(Array::Integer(Box::new(values))))
    };

    Some(value)
}

fn get_int32_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    use crate::record::codec::value::Int32;

    let src = src.get(range::<i32>(i, len))?;

    let value = if len == 1 {
        // SAFETY: `src` is 2 bytes.
        match Int32::from(i32::from_le_bytes(src.try_into().unwrap())) {
            Int32::Value(n) => Some(Value::Integer(n)),
            Int32::Missing => None,
            Int32::EndOfVector | Int32::Reserved(_) => todo!(),
        }
    } else {
        let values = Values::<'_, i32>::new(src);
        Some(Value::Array(Array::Integer(Box::new(values))))
    };

    Some(value)
}

fn get_float_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    use crate::record::codec::value::Float;

    let src = src.get(range::<f32>(i, len))?;

    let value = if len == 1 {
        // SAFETY: `src` is 2 bytes.
        match Float::from(f32::from_le_bytes(src.try_into().unwrap())) {
            Float::Value(n) => Some(Value::Float(n)),
            Float::Missing => None,
            Float::EndOfVector | Float::Reserved(_) => todo!(),
        }
    } else {
        let values = Values::<'_, f32>::new(src);
        Some(Value::Array(Array::Float(Box::new(values))))
    };

    Some(value)
}

fn get_string_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    const NUL: u8 = 0x00;

    let src = src.get(range::<u8>(i, len))?;

    let src = match src.iter().position(|&b| b == NUL) {
        Some(i) => &src[..i],
        None => src,
    };

    let s = str::from_utf8(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .unwrap(); // TODO

    Some(Some(Value::String(s)))
}
