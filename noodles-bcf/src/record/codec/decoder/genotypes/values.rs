use std::{error, fmt, mem, str};

use noodles_vcf::record::genotypes::sample::Value;

use crate::{
    lazy::record::value::{Float, Int16, Int32, Int8, Type},
    record::codec::decoder::value::{read_type, ty},
};

pub(super) fn read_values(
    src: &mut &[u8],
    sample_count: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    match read_type(src).map_err(DecodeError::InvalidType)? {
        Some(Type::Int8(0)) => Err(DecodeError::InvalidLength),
        Some(Type::Int8(1)) => read_int8_values(src, sample_count),
        Some(Type::Int8(n)) => read_int8_array_values(src, sample_count, n),
        Some(Type::Int16(0)) => Err(DecodeError::InvalidLength),
        Some(Type::Int16(1)) => read_int16_values(src, sample_count),
        Some(Type::Int16(n)) => read_int16_array_values(src, sample_count, n),
        Some(Type::Int32(0)) => Err(DecodeError::InvalidLength),
        Some(Type::Int32(1)) => read_int32_values(src, sample_count),
        Some(Type::Int32(n)) => read_int32_array_values(src, sample_count, n),
        Some(Type::Float(0)) => Err(DecodeError::InvalidLength),
        Some(Type::Float(1)) => read_float_values(src, sample_count),
        Some(Type::Float(n)) => read_float_array_values(src, sample_count, n),
        Some(Type::String(n)) => read_string_values(src, sample_count, n),
        ty => todo!("unhandled type: {ty:?}"),
    }
}

fn read_int8_values(
    src: &mut &[u8],
    sample_count: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = read_i8(src).map(Int8::from)?;

        match value {
            Int8::Value(n) => values.push(Some(Value::from(i32::from(n)))),
            Int8::Missing => values.push(None),
            _ => todo!("unhandled i8 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_int8_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let buf = read_i8s(src, len)?;

        let vs: Vec<_> = buf
            .into_iter()
            .map(Int8::from)
            .filter_map(|value| match value {
                Int8::Value(n) => Some(Some(i32::from(n))),
                Int8::Missing => Some(None),
                Int8::EndOfVector => None,
                _ => todo!("unhandled i8 array value: {:?}", value),
            })
            .collect();

        if vs.len() == 1 && vs[0].is_none() {
            values.push(None);
        } else {
            values.push(Some(Value::from(vs)));
        }
    }

    Ok(values)
}

fn read_int16_values(
    src: &mut &[u8],
    sample_count: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = read_i16(src).map(Int16::from)?;

        match value {
            Int16::Value(n) => values.push(Some(Value::from(i32::from(n)))),
            Int16::Missing => values.push(None),
            _ => todo!("unhandled i16 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_int16_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let buf = read_i16s(src, len)?;

        let vs: Vec<_> = buf
            .into_iter()
            .map(Int16::from)
            .filter_map(|value| match value {
                Int16::Value(n) => Some(Some(i32::from(n))),
                Int16::Missing => Some(None),
                Int16::EndOfVector => None,
                _ => todo!("unhandled i16 array value: {:?}", value),
            })
            .collect();

        if vs.len() == 1 && vs[0].is_none() {
            values.push(None);
        } else {
            values.push(Some(Value::from(vs)));
        }
    }

    Ok(values)
}

fn read_int32_values(
    src: &mut &[u8],
    sample_count: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = read_i32(src).map(Int32::from)?;

        match value {
            Int32::Value(n) => values.push(Some(Value::from(n))),
            Int32::Missing => values.push(None),
            _ => todo!("unhandled i32 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_int32_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let buf = read_i32s(src, len)?;

        let vs: Vec<_> = buf
            .into_iter()
            .map(Int32::from)
            .filter_map(|value| match value {
                Int32::Value(n) => Some(Some(n)),
                Int32::Missing => Some(None),
                Int32::EndOfVector => None,
                _ => todo!("unhandled i32 array value: {:?}", value),
            })
            .collect();

        if vs.len() == 1 && vs[0].is_none() {
            values.push(None);
        } else {
            values.push(Some(Value::from(vs)));
        }
    }

    Ok(values)
}

fn read_float_values(
    src: &mut &[u8],
    sample_count: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = read_f32(src).map(Float::from)?;

        match value {
            Float::Value(n) => values.push(Some(Value::from(n))),
            Float::Missing => values.push(None),
            _ => todo!("unhandled f32 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_float_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let buf = read_f32s(src, len)?;

        let vs: Vec<_> = buf
            .into_iter()
            .map(Float::from)
            .filter_map(|value| match value {
                Float::Value(n) => Some(Some(n)),
                Float::Missing => Some(None),
                Float::EndOfVector => None,
                _ => todo!("unhandled f32 array value: {:?}", value),
            })
            .collect();

        if vs.len() == 1 && vs[0].is_none() {
            values.push(None);
        } else {
            values.push(Some(Value::from(vs)));
        }
    }

    Ok(values)
}

fn read_string_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    const NUL: u8 = 0x00;

    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let buf = read_string(src, len)?;

        let data = match buf.iter().position(|&b| b == NUL) {
            Some(i) => &buf[..i],
            None => buf,
        };

        let s = str::from_utf8(data).map_err(DecodeError::InvalidString)?;
        let value = Value::from(s);

        values.push(Some(value));
    }

    Ok(values)
}

pub(super) fn read_genotype_values(
    src: &mut &[u8],
    sample_count: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    match read_type(src).map_err(DecodeError::InvalidType)? {
        Some(Type::Int8(len)) => match len {
            0 => values.push(None),
            1 => {
                for _ in 0..sample_count {
                    let value = read_i8(src)
                        .map(|v| parse_genotype_values(&[v]))
                        .map(Value::from)?;

                    values.push(Some(value));
                }
            }
            _ => {
                for _ in 0..sample_count {
                    let buf = read_i8s(src, len)?;
                    let value = Value::from(parse_genotype_values(&buf));
                    values.push(Some(value));
                }
            }
        },
        ty => todo!("unhandled type: {:?}", ty),
    }

    Ok(values)
}

fn parse_genotype_values(values: &[i8]) -> String {
    use std::fmt::Write;

    let mut genotype = String::new();

    for (i, &value) in values.iter().enumerate() {
        if let Int8::EndOfVector = Int8::from(value) {
            break;
        }

        let j = (value >> 1) - 1;
        let is_phased = value & 0x01 == 1;

        if i > 0 {
            if is_phased {
                genotype.push('|');
            } else {
                genotype.push('/');
            }
        }

        if j == -1 {
            genotype.push('.');
        } else {
            let _ = write!(genotype, "{j}");
        }
    }

    genotype
}

fn read_i8(src: &mut &[u8]) -> Result<i8, DecodeError> {
    if let Some((b, rest)) = src.split_first() {
        *src = rest;
        Ok(*b as i8)
    } else {
        Err(DecodeError::UnexpectedEof)
    }
}

fn read_i8s(src: &mut &[u8], len: usize) -> Result<Vec<i8>, DecodeError> {
    if src.len() < len {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(len);
    let values = buf.iter().map(|&b| b as i8).collect();
    *src = rest;

    Ok(values)
}

fn read_i16(src: &mut &[u8]) -> Result<i16, DecodeError> {
    if src.len() < mem::size_of::<i16>() {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(mem::size_of::<i16>());

    // SAFETY: `buf` is 2 bytes.
    let n = i16::from_le_bytes(buf.try_into().unwrap());

    *src = rest;

    Ok(n)
}

fn read_i16s(src: &mut &[u8], len: usize) -> Result<Vec<i16>, DecodeError> {
    let len = mem::size_of::<i16>() * len;

    if src.len() < len {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(len);

    let values = buf
        .chunks_exact(mem::size_of::<i16>())
        .map(|chunk| {
            // SAFETY: `chunk` is 2 bytes.
            i16::from_le_bytes(chunk.try_into().unwrap())
        })
        .collect();

    *src = rest;

    Ok(values)
}

fn read_i32(src: &mut &[u8]) -> Result<i32, DecodeError> {
    if src.len() < mem::size_of::<i32>() {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(mem::size_of::<i32>());

    // SAFETY: `buf` is 4 bytes.
    let n = i32::from_le_bytes(buf.try_into().unwrap());

    *src = rest;

    Ok(n)
}

fn read_i32s(src: &mut &[u8], len: usize) -> Result<Vec<i32>, DecodeError> {
    let len = mem::size_of::<i32>() * len;

    if src.len() < len {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(len);

    let values = buf
        .chunks_exact(mem::size_of::<i32>())
        .map(|chunk| {
            // SAFETY: `chunk` is 4 bytes.
            i32::from_le_bytes(chunk.try_into().unwrap())
        })
        .collect();

    *src = rest;

    Ok(values)
}

fn read_f32(src: &mut &[u8]) -> Result<f32, DecodeError> {
    if src.len() < mem::size_of::<f32>() {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(mem::size_of::<f32>());

    // SAFETY: `buf` is 4 bytes.
    let n = f32::from_le_bytes(buf.try_into().unwrap());

    *src = rest;

    Ok(n)
}

fn read_f32s(src: &mut &[u8], len: usize) -> Result<Vec<f32>, DecodeError> {
    let len = mem::size_of::<f32>() * len;

    if src.len() < len {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(len);

    let values = buf
        .chunks_exact(mem::size_of::<f32>())
        .map(|chunk| {
            // SAFETY: `chunk` is 4 bytes.
            f32::from_le_bytes(chunk.try_into().unwrap())
        })
        .collect();

    *src = rest;

    Ok(values)
}

fn read_string<'a>(src: &mut &'a [u8], len: usize) -> Result<&'a [u8], DecodeError> {
    if src.len() < len {
        return Err(DecodeError::UnexpectedEof);
    }

    let (buf, rest) = src.split_at(len);
    *src = rest;

    Ok(buf)
}

#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    UnexpectedEof,
    InvalidType(ty::DecodeError),
    InvalidLength,
    InvalidString(str::Utf8Error),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidType(e) => Some(e),
            Self::InvalidString(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidType(_) => write!(f, "invalid type"),
            Self::InvalidLength => write!(f, "invalid length"),
            Self::InvalidString(_) => write!(f, "invalid string"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_values_with_int8_values() {
        let mut src = &[
            0x11, // Some(Type::Int8(1))
            0x05, // Some(5)
            0x08, // Some(8)
            0x80, // None
        ][..];

        assert_eq!(
            read_values(&mut src, 3),
            Ok(vec![Some(Value::from(5)), Some(Value::from(8)), None])
        );
    }

    #[test]
    fn test_read_values_with_int8_array_values() {
        let mut src = &[
            0x21, // Some(Type::Int8(2))
            0x05, 0x08, // Some([Some(5), Some(8)])
            0x0d, 0x80, // Some([Some(13), None])
            0x15, 0x81, // Some([Some(21)])
            0x80, 0x81, // None
        ][..];

        assert_eq!(
            read_values(&mut src, 4),
            Ok(vec![
                Some(Value::from(vec![Some(5), Some(8)])),
                Some(Value::from(vec![Some(13), None])),
                Some(Value::from(vec![Some(21)])),
                None,
            ])
        );
    }

    #[test]
    fn test_read_values_with_int16_values() {
        let mut src = &[
            0x12, // Some(Type::Int16(1))
            0x05, 0x00, // Some(5)
            0x08, 0x00, // Some(8)
            0x00, 0x80, // None
        ][..];

        assert_eq!(
            read_values(&mut src, 3),
            Ok(vec![Some(Value::from(5)), Some(Value::from(8)), None])
        );
    }

    #[test]
    fn test_read_values_with_int16_array_values() {
        let mut src = &[
            0x22, // Some(Type::Int16(2))
            0x05, 0x00, 0x08, 0x00, // Some([Some(5), Some(8)])
            0x0d, 0x00, 0x00, 0x80, // Some([Some(13), None])
            0x15, 0x00, 0x01, 0x80, // Some([Some(21)])
            0x00, 0x80, 0x01, 0x80, // None
        ][..];

        assert_eq!(
            read_values(&mut src, 4),
            Ok(vec![
                Some(Value::from(vec![Some(5), Some(8)])),
                Some(Value::from(vec![Some(13), None])),
                Some(Value::from(vec![Some(21)])),
                None,
            ])
        );
    }

    #[test]
    fn test_read_values_with_int32_values() {
        let mut src = &[
            0x13, // Some(Type::Int32(1))
            0x05, 0x00, 0x00, 0x00, // Some(5)
            0x08, 0x00, 0x00, 0x00, // Some(8)
            0x00, 0x00, 0x00, 0x80, // None
        ][..];

        assert_eq!(
            read_values(&mut src, 3),
            Ok(vec![Some(Value::from(5)), Some(Value::from(8)), None])
        );
    }

    #[test]
    fn test_read_values_with_int32_array_values() {
        let mut src = &[
            0x23, // Some(Type::Int32(2))
            0x05, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x00, // Some([Some(5), Some(8)])
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, // Some([Some(13), None])
            0x15, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, // Some([Some(21)])
            0x00, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // None
        ][..];

        assert_eq!(
            read_values(&mut src, 4),
            Ok(vec![
                Some(Value::from(vec![Some(5), Some(8)])),
                Some(Value::from(vec![Some(13), None])),
                Some(Value::from(vec![Some(21)])),
                None,
            ])
        );
    }

    #[test]
    fn test_read_values_with_float_values() {
        let mut src = &[
            0x15, // Some(Type::Float(1))
            0x00, 0x00, 0x00, 0x00, // Some(0.0)
            0x00, 0x00, 0x80, 0x3f, // Some(1.0)
            0x01, 0x00, 0x80, 0x7f, // None
        ][..];

        assert_eq!(
            read_values(&mut src, 3),
            Ok(vec![Some(Value::from(0.0)), Some(Value::from(1.0)), None])
        );
    }

    #[test]
    fn test_read_values_with_float_array_values() {
        let mut src = &[
            0x25, // Some(Type::Float(2))
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x3f, // Some([Some(0.0), Some(1.0)])
            0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x80, 0x7f, // Some([Some(0.0), None])
            0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x80, 0x7f, // Some([Some(0.0)])
            0x01, 0x00, 0x80, 0x7f, 0x02, 0x00, 0x80, 0x7f, // None
        ][..];

        assert_eq!(
            read_values(&mut src, 4),
            Ok(vec![
                Some(Value::from(vec![Some(0.0), Some(1.0)])),
                Some(Value::from(vec![Some(0.0), None])),
                Some(Value::from(vec![Some(0.0)])),
                None,
            ])
        );
    }

    #[test]
    fn test_read_values_with_string_values() {
        let mut src = &[
            0x47, // Some(Type::String(4))
            b'n', 0x00, 0x00, 0x00, // "n"
            b'n', b'd', b'l', 0x00, // "ndl"
            b'n', b'd', b'l', b's', // "ndls"
        ][..];

        assert_eq!(
            read_values(&mut src, 3),
            Ok(vec![
                Some(Value::from("n")),
                Some(Value::from("ndl")),
                Some(Value::from("ndls")),
            ])
        );
    }

    #[test]
    fn test_parse_genotype_genotype_field_values() {
        // Examples from § 6.3.3 Type encoding (2021-05-13)

        assert_eq!(parse_genotype_values(&[0x02, 0x02]), "0/0");
        assert_eq!(parse_genotype_values(&[0x02, 0x04]), "0/1");
        assert_eq!(parse_genotype_values(&[0x04, 0x04]), "1/1");
        assert_eq!(parse_genotype_values(&[0x02, 0x05]), "0|1");
        assert_eq!(parse_genotype_values(&[0x00, 0x00]), "./.");
        assert_eq!(parse_genotype_values(&[0x02]), "0");
        assert_eq!(parse_genotype_values(&[0x04]), "1");
        assert_eq!(parse_genotype_values(&[0x02, 0x04, 0x06]), "0/1/2");
        assert_eq!(parse_genotype_values(&[0x02, 0x04, 0x07]), "0/1|2");
        assert_eq!(
            parse_genotype_values(&[0x02, i8::from(Int8::EndOfVector)]),
            "0"
        );
    }
}
