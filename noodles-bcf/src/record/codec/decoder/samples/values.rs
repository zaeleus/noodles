use std::{error, fmt, str};

use memchr::memchr;
use noodles_vcf::{
    header::record::value::map::format::{self, Number},
    variant::record_buf::samples::sample::{Value, value::Genotype},
};

use crate::record::codec::{
    decoder::{
        raw_value::{
            self, read_f32, read_f32s, read_i8, read_i8s, read_i16, read_i16s, read_i32, read_i32s,
            read_string,
        },
        value::{read_type, ty},
    },
    value::{Float, Int8, Int16, Int32, Type},
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
        (_, format::Type::Character | format::Type::String, Type::String(0)) => {
            Ok(vec![None; sample_count])
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
        (_, format::Type::Character, Type::String(n)) => {
            read_char_array_values(src, sample_count, n)
        }
        (_, format::Type::String, Type::String(n)) => {
            read_string_array_values(src, sample_count, n)
        }

        _ => Err(DecodeError::TypeMismatch(value_ty, ty)),
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
            Int8::EndOfVector | Int8::Reserved(_) => return Err(DecodeError::InvalidInteger),
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
            .map_while(|value| match Int8::from(value) {
                Int8::Value(n) => Some(Ok(Some(i32::from(n)))),
                Int8::Missing => Some(Ok(None)),
                Int8::EndOfVector => None,
                Int8::Reserved(_) => Some(Err(DecodeError::InvalidInteger)),
            })
            .collect::<Result<_, _>>()?;

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
            Int16::EndOfVector | Int16::Reserved(_) => return Err(DecodeError::InvalidInteger),
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
            .map_while(|value| match Int16::from(value) {
                Int16::Value(n) => Some(Ok(Some(i32::from(n)))),
                Int16::Missing => Some(Ok(None)),
                Int16::EndOfVector => None,
                Int16::Reserved(_) => Some(Err(DecodeError::InvalidInteger)),
            })
            .collect::<Result<_, _>>()?;

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
            Int32::EndOfVector | Int32::Reserved(_) => return Err(DecodeError::InvalidInteger),
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
            .map_while(|value| match Int32::from(value) {
                Int32::Value(n) => Some(Ok(Some(n))),
                Int32::Missing => Some(Ok(None)),
                Int32::EndOfVector => None,
                Int32::Reserved(_) => Some(Err(DecodeError::InvalidInteger)),
            })
            .collect::<Result<_, _>>()?;

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
            Float::EndOfVector | Float::Reserved(_) => return Err(DecodeError::InvalidFloat),
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
            .map_while(|value| match Float::from(value) {
                Float::Value(n) => Some(Ok(Some(n))),
                Float::Missing => Some(Ok(None)),
                Float::EndOfVector => None,
                Float::Reserved(_) => Some(Err(DecodeError::InvalidFloat)),
            })
            .collect::<Result<_, _>>()?;

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

    let data = match memchr(NUL, buf) {
        Some(i) => &buf[..i],
        None => buf,
    };

    str::from_utf8(data).map_err(DecodeError::InvalidString)
}

fn read_nonempty_string_until_nul<'a>(
    src: &mut &'a [u8],
    len: usize,
) -> Result<&'a str, DecodeError> {
    let s = read_string_until_nul(src, len)?;

    if s.is_empty() {
        Err(DecodeError::EmptyValue)
    } else {
        Ok(s)
    }
}

fn read_char_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    const MISSING: char = '.';

    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let s = read_nonempty_string_until_nul(src, len)?;
        let mut chars = s.chars();

        let value = match (chars.next(), chars.next()) {
            (Some(MISSING), None) => None,
            (Some(c), None) => Some(Value::from(c)),
            _ => return Err(DecodeError::InvalidCharacter),
        };

        values.push(value);
    }

    Ok(values)
}

fn read_char_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    const DELIMITER: char = ',';
    const MISSING: &str = ".";
    const MISSING_VALUE: char = '.';

    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = match read_nonempty_string_until_nul(src, len)? {
            MISSING => None,
            s => {
                let cs: Vec<_> = s
                    .split(DELIMITER)
                    .map(|t| {
                        if t.is_empty() {
                            return Err(DecodeError::EmptyValue);
                        }

                        let mut chars = t.chars();

                        match (chars.next(), chars.next()) {
                            (Some(MISSING_VALUE), None) => Ok(None),
                            (Some(c), None) => Ok(Some(c)),
                            _ => Err(DecodeError::InvalidCharacter),
                        }
                    })
                    .collect::<Result<_, _>>()?;

                Some(Value::from(cs))
            }
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
    const MISSING: &str = ".";

    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = match read_nonempty_string_until_nul(src, len)? {
            MISSING => None,
            s => Some(Value::from(s)),
        };

        values.push(value);
    }

    Ok(values)
}

fn read_string_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    const DELIMITER: char = ',';
    const MISSING: &str = ".";
    const EMPTY: &str = "";

    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = match read_nonempty_string_until_nul(src, len)? {
            MISSING => None,
            s => {
                let vs: Vec<Option<String>> = s
                    .split(DELIMITER)
                    .map(|t| match t {
                        EMPTY => Err(DecodeError::EmptyValue),
                        MISSING => Ok(None),
                        _ => Ok(Some(t.into())),
                    })
                    .collect::<Result<_, _>>()?;

                Some(Value::from(vs))
            }
        };

        values.push(value);
    }

    Ok(values)
}

enum GenotypeValues<'a> {
    Int8(&'a [i8]),
    Int16(&'a [i16]),
    Int32(&'a [i32]),
}

impl GenotypeValues<'_> {
    fn iter(&self) -> Box<dyn Iterator<Item = Result<i32, DecodeError>> + '_> {
        match self {
            Self::Int8(values) => Box::new(values.iter().map_while(|n| match Int8::from(*n) {
                Int8::Value(n) => Some(Ok(i32::from(n))),
                Int8::EndOfVector => None,
                Int8::Missing | Int8::Reserved(_) => {
                    Some(Err(DecodeError::InvalidGenotypeSampleValue))
                }
            })),
            Self::Int16(values) => Box::new(values.iter().map_while(|n| match Int16::from(*n) {
                Int16::Value(n) => Some(Ok(i32::from(n))),
                Int16::EndOfVector => None,
                Int16::Missing | Int16::Reserved(_) => {
                    Some(Err(DecodeError::InvalidGenotypeSampleValue))
                }
            })),
            Self::Int32(values) => Box::new(values.iter().map_while(|n| match Int32::from(*n) {
                Int32::Value(n) => Some(Ok(n)),
                Int32::EndOfVector => None,
                Int32::Missing | Int32::Reserved(_) => {
                    Some(Err(DecodeError::InvalidGenotypeSampleValue))
                }
            })),
        }
    }
}

pub(super) fn read_genotype_values(
    src: &mut &[u8],
    sample_count: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    match read_type(src).map_err(DecodeError::InvalidType)? {
        Some(Type::Int8(len)) => read_i8_genotype_values(src, sample_count, len),
        Some(Type::Int16(len)) => read_i16_genotype_values(src, sample_count, len),
        Some(Type::Int32(len)) => read_i32_genotype_values(src, sample_count, len),
        _ => Err(DecodeError::InvalidGenotypeType),
    }
}

fn read_i8_genotype_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    match len {
        0 => values.push(None),
        1 => {
            for _ in 0..sample_count {
                let n = read_i8(src).map_err(DecodeError::InvalidRawValue)?;
                let vs = GenotypeValues::Int8(&[n]);
                let genotype = parse_genotype_values(vs)?;
                values.push(Some(Value::Genotype(genotype)));
            }
        }
        _ => {
            for _ in 0..sample_count {
                let buf = read_i8s(src, len).map_err(DecodeError::InvalidRawValue)?;
                let vs = GenotypeValues::Int8(&buf);
                let genotype = parse_genotype_values(vs)?;
                values.push(Some(Value::Genotype(genotype)));
            }
        }
    }

    Ok(values)
}

fn read_i16_genotype_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    match len {
        0 => values.push(None),
        1 => {
            for _ in 0..sample_count {
                let n = read_i16(src).map_err(DecodeError::InvalidRawValue)?;
                let vs = GenotypeValues::Int16(&[n]);
                let genotype = parse_genotype_values(vs)?;
                values.push(Some(Value::Genotype(genotype)));
            }
        }
        _ => {
            for _ in 0..sample_count {
                let buf = read_i16s(src, len).map_err(DecodeError::InvalidRawValue)?;
                let vs = GenotypeValues::Int16(&buf);
                let genotype = parse_genotype_values(vs)?;
                values.push(Some(Value::Genotype(genotype)));
            }
        }
    }

    Ok(values)
}

fn read_i32_genotype_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> Result<Vec<Option<Value>>, DecodeError> {
    let mut values = Vec::with_capacity(sample_count);

    match len {
        0 => values.push(None),
        1 => {
            for _ in 0..sample_count {
                let n = read_i32(src).map_err(DecodeError::InvalidRawValue)?;
                let vs = GenotypeValues::Int32(&[n]);
                let genotype = parse_genotype_values(vs)?;
                values.push(Some(Value::Genotype(genotype)));
            }
        }
        _ => {
            for _ in 0..sample_count {
                let buf = read_i32s(src, len).map_err(DecodeError::InvalidRawValue)?;
                let vs = GenotypeValues::Int32(&buf);
                let genotype = parse_genotype_values(vs)?;
                values.push(Some(Value::Genotype(genotype)));
            }
        }
    }

    Ok(values)
}

fn parse_genotype_values(values: GenotypeValues) -> Result<Genotype, DecodeError> {
    use noodles_vcf::variant::{
        record::samples::series::value::genotype::Phasing,
        record_buf::samples::sample::value::genotype::Allele,
    };

    let mut alleles = Vec::new();

    for result in values.iter() {
        let value = result?;

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
    TypeMismatch(Type, format::Type),
    InvalidRawValue(raw_value::DecodeError),
    InvalidInteger,
    InvalidFloat,
    InvalidCharacter,
    InvalidString(str::Utf8Error),
    EmptyValue,
    InvalidGenotype,
    InvalidGenotypeType,
    InvalidGenotypeSampleValue,
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
            Self::TypeMismatch(actual, expected) => {
                let resolved_actual = match actual {
                    Type::Int8(_) | Type::Int16(_) | Type::Int32(_) => format::Type::Integer,
                    Type::Float(_) => format::Type::Float,
                    Type::String(_) => format::Type::String,
                };

                write!(
                    f,
                    "type mismatch: expected {expected:?}, got {resolved_actual:?}"
                )
            }
            Self::InvalidRawValue(_) => write!(f, "invalid raw value"),
            Self::InvalidInteger => write!(f, "invalid integer"),
            Self::InvalidFloat => write!(f, "invalid float"),
            Self::InvalidCharacter => write!(f, "invalid character"),
            Self::InvalidString(_) => write!(f, "invalid string"),
            Self::EmptyValue => write!(f, "empty value"),
            Self::InvalidGenotype => write!(f, "invalid genotype"),
            Self::InvalidGenotypeType => write!(f, "invalid genotype type"),
            Self::InvalidGenotypeSampleValue => write!(f, "invalid genotype sample value"),
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

        let mut src = &[
            0x11, // Some(Type::Int8(1))
            0x81, // Int8::EndOfVector
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::Integer, 1),
            Err(DecodeError::InvalidInteger)
        );

        let mut src = &[
            0x11, // Some(Type::Int8(1))
            0x82, // Int8::Reserved(0x82)
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::Integer, 1),
            Err(DecodeError::InvalidInteger)
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

        let mut src = &[
            0x21, // Some(Type::Int8(2))
            0x82, 0x81, // [Int8::Reserved(0x82)]
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::Integer, 1),
            Err(DecodeError::InvalidInteger)
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

        let mut src = &[
            0x12, // Some(Type::Int16(1))
            0x01, 0x80, // Int16::EndOfVector
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::Integer, 1),
            Err(DecodeError::InvalidInteger)
        );

        let mut src = &[
            0x12, // Some(Type::Int16(1))
            0x02, 0x80, // Int16::Reserved(0x8002)
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::Integer, 1),
            Err(DecodeError::InvalidInteger)
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

        let mut src = &[
            0x22, // Some(Type::Int16(2))
            0x02, 0x80, 0x01, 0x80, // Int16::Reserved(0x8002)]
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::Integer, 1),
            Err(DecodeError::InvalidInteger)
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

        let mut src = &[
            0x13, // Some(Type::Int32(1))
            0x01, 0x00, 0x00, 0x80, // Int32::EndOfVector
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::Integer, 1),
            Err(DecodeError::InvalidInteger)
        );

        let mut src = &[
            0x13, // Some(Type::Int32(1))
            0x02, 0x00, 0x00, 0x80, // Int32::Reserved(0x80000002)
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::Integer, 1),
            Err(DecodeError::InvalidInteger)
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

        let mut src = &[
            0x23, // Some(Type::Int32(2))
            0x02, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // [Int32::Reserved(0x80000002)]
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::Integer, 1),
            Err(DecodeError::InvalidInteger)
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

        let mut src = &[
            0x15, // Some(Type::Float(1))
            0x02, 0x00, 0x80, 0x7f, // Float::EndOfVector
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::Float, 1),
            Err(DecodeError::InvalidFloat)
        );

        let mut src = &[
            0x15, // Some(Type::Float(1))
            0x03, 0x00, 0x80, 0x7f, // Float::Reserved(0x7f800003)
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::Float, 1),
            Err(DecodeError::InvalidFloat)
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

        let mut src = &[
            0x25, // Some(Type::Float(2))
            0x03, 0x00, 0x80, 0x7f, 0x02, 0x00, 0x80, 0x7f, // [Float::Reserved(0x7f800003)]
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::Float, 1),
            Err(DecodeError::InvalidFloat)
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

        let mut src = &[0x07][..]; // Some(Type::String(0))
        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::Character, 2),
            Ok(vec![None, None])
        );

        let mut src = &[
            0x17, // Some(Type::String(1))
            0x00, // NUL
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::Character, 1),
            Err(DecodeError::EmptyValue)
        );

        let mut src = &[
            0x27, // Some(Type::String(2))
            b'n', b'd', // "nd"
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::Character, 1),
            Err(DecodeError::InvalidCharacter)
        );
    }

    #[test]
    fn test_read_values_with_character_array_values() {
        let mut src = &[
            0x37, // Some(Type::String(3))
            b'n', b',', b'd', // [Some('n'), Some('d')]
            b'n', b',', b'.', // [Some('n'), None]
            b'n', 0x00, 0x00, // [Some('n')]
            b'.', 0x00, 0x00, // None
        ][..];

        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::Character, 4),
            Ok(vec![
                Some(Value::from(vec![Some('n'), Some('d')])),
                Some(Value::from(vec![Some('n'), None])),
                Some(Value::from(vec![Some('n')])),
                None,
            ])
        );

        let mut src = &[0x07][..]; // Some(Type::String(0))
        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::Character, 2),
            Ok(vec![None, None])
        );

        let mut src = &[
            0x37, // Some(Type::String(3))
            0x00, 0x00, 0x00, // NUL
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::Character, 1),
            Err(DecodeError::EmptyValue)
        );

        let mut src = &[
            0x37, // Some(Type::String(3))
            b'n', b',', // "n"
            0x00, // NUL
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::Character, 1),
            Err(DecodeError::EmptyValue)
        );

        let mut src = &[
            0x47, // Some(Type::String(4))
            b'n', b',', // "n"
            b'd', b'l', // "dl"
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::Character, 1),
            Err(DecodeError::InvalidCharacter)
        );
    }

    #[test]
    fn test_read_values_with_string_values() {
        let mut src = &[
            0x47, // Some(Type::String(4))
            b'n', b'd', b'l', b's', // Some("ndls")
            b'n', 0x00, 0x00, 0x00, // Some("n")
            b'.', 0x00, 0x00, 0x00, // None
        ][..];

        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::String, 3),
            Ok(vec![
                Some(Value::from("ndls")),
                Some(Value::from("n")),
                None,
            ])
        );

        let mut src = &[0x07][..]; // Some(Type::String(0))
        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::String, 2),
            Ok(vec![None, None])
        );

        let mut src = &[
            0x17, // Some(Type::String(1))
            0x00, // NUL
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(1), format::Type::String, 1),
            Err(DecodeError::EmptyValue)
        );
    }

    #[test]
    fn test_read_values_with_string_array_values() {
        let mut src = &[
            0x37, // Some(Type::String(3))
            b'n', b',', b'd', // [Some("n"), Some("d")]
            b'n', b',', b'.', // [Some("n"), None]
            b'n', 0x00, 0x00, // [Some("n")]
            b'.', 0x00, 0x00, // None
        ][..];

        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::String, 4),
            Ok(vec![
                Some(Value::from(vec![
                    Some(String::from("n")),
                    Some(String::from("d"))
                ])),
                Some(Value::from(vec![Some(String::from("n")), None])),
                Some(Value::from(vec![Some(String::from("n"))])),
                None,
            ])
        );

        let mut src = &[0x07][..]; // Some(Type::String(0))
        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::String, 2),
            Ok(vec![None, None])
        );

        let mut src = &[
            0x37, // Some(Type::String(3))
            0x00, 0x00, 0x00, // NUL
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::String, 1),
            Err(DecodeError::EmptyValue)
        );

        let mut src = &[
            0x37, // Some(Type::String(3))
            b'n', b',', // "n"
            0x00, // NUL
        ][..];
        assert_eq!(
            read_values(&mut src, Number::Count(2), format::Type::String, 1),
            Err(DecodeError::EmptyValue)
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
            parse_genotype_values(GenotypeValues::Int8(&[0x02, 0x02]))?,
            [
                Allele::new(Some(0), Phasing::Unphased),
                Allele::new(Some(0), Phasing::Unphased)
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(GenotypeValues::Int8(&[0x02, 0x04]))?,
            [
                Allele::new(Some(0), Phasing::Unphased),
                Allele::new(Some(1), Phasing::Unphased)
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(GenotypeValues::Int8(&[0x04, 0x04]))?,
            [
                Allele::new(Some(1), Phasing::Unphased),
                Allele::new(Some(1), Phasing::Unphased)
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(GenotypeValues::Int8(&[0x02, 0x05]))?,
            [
                Allele::new(Some(0), Phasing::Unphased),
                Allele::new(Some(1), Phasing::Phased)
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(GenotypeValues::Int8(&[0x00, 0x00]))?,
            [
                Allele::new(None, Phasing::Unphased),
                Allele::new(None, Phasing::Unphased)
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(GenotypeValues::Int8(&[0x02]))?,
            [Allele::new(Some(0), Phasing::Unphased)]
                .into_iter()
                .collect()
        );

        assert_eq!(
            parse_genotype_values(GenotypeValues::Int8(&[0x04]))?,
            [Allele::new(Some(1), Phasing::Unphased)]
                .into_iter()
                .collect()
        );

        assert_eq!(
            parse_genotype_values(GenotypeValues::Int8(&[0x02, 0x04, 0x06]))?,
            [
                Allele::new(Some(0), Phasing::Unphased),
                Allele::new(Some(1), Phasing::Unphased),
                Allele::new(Some(2), Phasing::Unphased),
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(GenotypeValues::Int8(&[0x02, 0x04, 0x07]))?,
            [
                Allele::new(Some(0), Phasing::Unphased),
                Allele::new(Some(1), Phasing::Unphased),
                Allele::new(Some(2), Phasing::Phased),
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            parse_genotype_values(GenotypeValues::Int8(&[0x02, i8::from(Int8::EndOfVector)]))?,
            [Allele::new(Some(0), Phasing::Unphased)]
                .into_iter()
                .collect()
        );

        Ok(())
    }
}
