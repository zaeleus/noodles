use std::{error, fmt, str};

use noodles_vcf::{
    header::record::value::map::format::{self, Number},
    variant::record_buf::samples::sample::{value::Genotype, Value},
};

use crate::record::codec::{
    decoder::{
        raw_value::{
            self, read_f32, read_f32s, read_i16, read_i16s, read_i32, read_i32s, read_i8, read_i8s,
            read_string,
        },
        value::{read_type, ty},
    },
    value::{Float, Int16, Int32, Int8, Type},
};

pub(super) fn read_values(
    src: &mut &[u8],
    number: Number,
    ty: format::Type,
    sample_count: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let value_ty = read_type(src)
        .map_err(DecodeError::InvalidType)?
        .expect("unhandled type");

    match (number, ty, value_ty) {
        (Number::Count(0), _, _) => todo!("invalid number for type"),

        (_, _, Type::Int8(0) | Type::Int16(0) | Type::Int32(0) | Type::Float(0)) => {
            Err(DecodeError::InvalidLength)
        }

        (Number::Count(1), format::Type::Integer, Type::Int8(1)) => {
            read_i8_values(src, sample_count)
        }
        (Number::Count(1), format::Type::Integer, Type::Int16(1)) => {
            read_i16_values(src, sample_count)
        }
        (Number::Count(1), format::Type::Integer, Type::Int32(1)) => {
            read_i32_values(src, sample_count)
        }
        (Number::Count(1), format::Type::Float, Type::Float(1)) => {
            read_f32_values(src, sample_count)
        }
        (Number::Count(1), format::Type::Character, Type::String(n)) => {
            read_char_values(src, sample_count, n)
        }
        (Number::Count(1), format::Type::String, Type::String(n)) => {
            read_string_values(src, sample_count, n)
        }

        (_, format::Type::Integer, Type::Int8(n)) => read_i8_array_values(src, sample_count, n),
        (_, format::Type::Integer, Type::Int16(n)) => read_i16_array_values(src, sample_count, n),
        (_, format::Type::Integer, Type::Int32(n)) => read_i32_array_values(src, sample_count, n),
        (_, format::Type::Float, Type::Float(n)) => read_f32_array_values(src, sample_count, n),
        (_, format::Type::Character, Type::String(_)) => todo!(),
        (_, format::Type::String, Type::String(n)) => {
            read_string_array_values(src, sample_count, n)
        }

        _ => todo!("unhandled type"),
    }
}

fn read_i8_values(src: &mut &[u8], sample_count: usize) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = read_i8(src)
            .map(Int8::from)
            .map_err(DecodeError::InvalidRawValue)?;

        match value {
            Int8::Value(n) => values.push(Some(Value::from(i32::from(n)))),
            Int8::Missing => values.push(None),
            _ => todo!("unhandled i8 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_i8_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let buf = read_i8s(src, len).map_err(DecodeError::InvalidRawValue)?;

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

fn read_i16_values(
    src: &mut &[u8],
    sample_count: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = read_i16(src)
            .map(Int16::from)
            .map_err(DecodeError::InvalidRawValue)?;

        match value {
            Int16::Value(n) => values.push(Some(Value::from(i32::from(n)))),
            Int16::Missing => values.push(None),
            _ => todo!("unhandled i16 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_i16_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let buf = read_i16s(src, len).map_err(DecodeError::InvalidRawValue)?;

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

fn read_i32_values(
    src: &mut &[u8],
    sample_count: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = read_i32(src)
            .map(Int32::from)
            .map_err(DecodeError::InvalidRawValue)?;

        match value {
            Int32::Value(n) => values.push(Some(Value::from(n))),
            Int32::Missing => values.push(None),
            _ => todo!("unhandled i32 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_i32_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let buf = read_i32s(src, len).map_err(DecodeError::InvalidRawValue)?;

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

fn read_f32_values(
    src: &mut &[u8],
    sample_count: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = read_f32(src)
            .map(Float::from)
            .map_err(DecodeError::InvalidRawValue)?;

        match value {
            Float::Value(n) => values.push(Some(Value::from(n))),
            Float::Missing => values.push(None),
            _ => todo!("unhandled f32 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_f32_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let buf = read_f32s(src, len).map_err(DecodeError::InvalidRawValue)?;

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

fn read_string_until_nul<'a>(src: &mut &'a [u8], len: usize) -> Result<&'a str, DecodeError> {
    const NUL: u8 = 0x00;

    let buf = read_string(src, len).map_err(DecodeError::InvalidRawValue)?;

    let data = match buf.iter().position(|&b| b == NUL) {
        Some(i) => &buf[..i],
        None => buf,
    };

    str::from_utf8(data).map_err(DecodeError::InvalidString)
}

fn read_char_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    const MISSING: char = '.';

    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let s = read_string_until_nul(src, len)?;
        let c = s.chars().next().unwrap();

        let value = match c {
            MISSING => None,
            _ => Some(Value::from(c)),
        };

        values.push(value);
    }

    Ok(values)
}

fn read_string_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let s = read_string_until_nul(src, len)?;
        let value = Value::from(s);
        values.push(Some(value));
    }

    Ok(values)
}

fn read_string_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    const NUL: u8 = 0x00;
    const DELIMITER: char = ',';

    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let buf = read_string(src, len).map_err(DecodeError::InvalidRawValue)?;

        let data = match buf.iter().position(|&b| b == NUL) {
            Some(i) => &buf[..i],
            None => buf,
        };

        let s = str::from_utf8(data).map_err(DecodeError::InvalidString)?;
        let ss: Vec<Option<String>> = s.split(DELIMITER).map(String::from).map(Some).collect();
        let value = Value::from(ss);

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
                        .map_err(DecodeError::InvalidRawValue)
                        .and_then(|v| parse_genotype_values(&[v]))
                        .map(Value::Genotype)?;

                    values.push(Some(value));
                }
            }
            _ => {
                for _ in 0..sample_count {
                    let buf = read_i8s(src, len).map_err(DecodeError::InvalidRawValue)?;
                    let genotype = parse_genotype_values(&buf)?;
                    let value = Value::Genotype(genotype);
                    values.push(Some(value));
                }
            }
        },
        ty => todo!("unhandled type: {:?}", ty),
    }

    Ok(values)
}

fn parse_genotype_values(values: &[i8]) -> Result<Genotype, DecodeError> {
    use noodles_vcf::variant::{
        record::samples::series::value::genotype::Phasing,
        record_buf::samples::sample::value::genotype::Allele,
    };

    let mut alleles = Vec::with_capacity(values.len());

    for &value in values {
        if let Int8::EndOfVector = Int8::from(value) {
            break;
        }

        let j = (value >> 1) - 1;
        let is_phased = value & 0x01 == 1;

        let phasing = if is_phased {
            Phasing::Phased
        } else {
            Phasing::Unphased
        };

        let position = if j == -1 {
            None
        } else {
            usize::try_from(j)
                .map(Some)
                .map_err(|_| DecodeError::InvalidGenotype)?
        };

        alleles.push(Allele::new(position, phasing));
    }

    Ok(alleles.into_iter().collect())
}

#[allow(clippy::enum_variant_names)]
#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    InvalidType(ty::DecodeError),
    InvalidLength,
    InvalidRawValue(raw_value::DecodeError),
    InvalidString(str::Utf8Error),
    InvalidGenotype,
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidType(e) => Some(e),
            Self::InvalidRawValue(e) => Some(e),
            Self::InvalidString(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidType(_) => write!(f, "invalid type"),
            Self::InvalidLength => write!(f, "invalid length"),
            Self::InvalidRawValue(_) => write!(f, "invalid raw value"),
            Self::InvalidString(_) => write!(f, "invalid string"),
            Self::InvalidGenotype => write!(f, "invalid genotype"),
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
            read_values(&mut src, Number::Count(1), format::Type::Integer, 3),
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
            read_values(&mut src, Number::Count(2), format::Type::Integer, 4),
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
            read_values(&mut src, Number::Count(1), format::Type::Integer, 3),
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
            read_values(&mut src, Number::Count(2), format::Type::Integer, 4),
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
            read_values(&mut src, Number::Count(1), format::Type::Integer, 3),
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
            read_values(&mut src, Number::Count(2), format::Type::Integer, 4),
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
            read_values(&mut src, Number::Count(1), format::Type::Float, 3),
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
            read_values(&mut src, Number::Count(2), format::Type::Float, 4),
            Ok(vec![
                Some(Value::from(vec![Some(0.0), Some(1.0)])),
                Some(Value::from(vec![Some(0.0), None])),
                Some(Value::from(vec![Some(0.0)])),
                None,
            ])
        );
    }

    #[test]
    fn test_read_values_with_character_values() {
        let mut src = &[
            0x17, // Some(Type::String(1))
            b'n', // Some('n')
            b'.', // None
        ][..];

        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::Character, 2),
            Ok(vec![Some(Value::from('n')), None])
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
            read_values(&mut src, Number::Count(1), format::Type::String, 3),
            Ok(vec![
                Some(Value::from("n")),
                Some(Value::from("ndl")),
                Some(Value::from("ndls")),
            ])
        );
    }

    #[test]
    fn test_read_values_with_string_array_values() {
        let mut src = &[
            0x47, // Some(Type::String(4))
            b'n', 0x00, 0x00, 0x00, // "n"
            b'n', b',', b'l', 0x00, // "n,l"
            b'n', b',', b'l', b's', // "n,ls"
        ][..];

        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::String, 3),
            Ok(vec![
                Some(Value::from(vec![Some(String::from("n"))])),
                Some(Value::from(vec![
                    Some(String::from("n")),
                    Some(String::from("l"))
                ])),
                Some(Value::from(vec![
                    Some(String::from("n")),
                    Some(String::from("ls"))
                ])),
            ])
        );
    }

    #[test]
    fn test_parse_genotype_values() -> Result<(), DecodeError> {
        use noodles_vcf::variant::{
            record::samples::series::value::genotype::Phasing,
            record_buf::samples::sample::value::genotype::Allele,
        };

        // Examples from § 6.3.3 "Type encoding" (2024-04-20).

        assert_eq!(
            parse_genotype_values(&[0x02, 0x02])?,
            [
                Allele::new(Some(0), Phasing::Unphased),
                Allele::new(Some(0), Phasing::Unphased)
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(&[0x02, 0x04])?,
            [
                Allele::new(Some(0), Phasing::Unphased),
                Allele::new(Some(1), Phasing::Unphased)
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(&[0x04, 0x04])?,
            [
                Allele::new(Some(1), Phasing::Unphased),
                Allele::new(Some(1), Phasing::Unphased)
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(&[0x02, 0x05])?,
            [
                Allele::new(Some(0), Phasing::Unphased),
                Allele::new(Some(1), Phasing::Phased)
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(&[0x00, 0x00])?,
            [
                Allele::new(None, Phasing::Unphased),
                Allele::new(None, Phasing::Unphased)
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(&[0x02])?,
            [Allele::new(Some(0), Phasing::Unphased)]
                .into_iter()
                .collect()
        );

        assert_eq!(
            parse_genotype_values(&[0x04])?,
            [Allele::new(Some(1), Phasing::Unphased)]
                .into_iter()
                .collect()
        );

        assert_eq!(
            parse_genotype_values(&[0x02, 0x04, 0x06])?,
            [
                Allele::new(Some(0), Phasing::Unphased),
                Allele::new(Some(1), Phasing::Unphased),
                Allele::new(Some(2), Phasing::Unphased),
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(&[0x02, 0x04, 0x07])?,
            [
                Allele::new(Some(0), Phasing::Unphased),
                Allele::new(Some(1), Phasing::Unphased),
                Allele::new(Some(2), Phasing::Phased),
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(&[0x02, i8::from(Int8::EndOfVector)])?,
            [Allele::new(Some(0), Phasing::Unphased)]
                .into_iter()
                .collect()
        );

        Ok(())
    }
}
