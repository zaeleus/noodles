use std::{
    borrow::Cow,
    cmp,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_vcf::{
    header::record::value::{Map, map::Format},
    variant::record::samples::series::{
        Value,
        value::Array,
        value::{Genotype, genotype::Phasing},
    },
};

use crate::{
    io::writer::num::write_i32_le,
    record::codec::{
        encoder::value::write_type,
        value::{Float, Int8, Int16, Int32, Type},
    },
};

const DELIMITER: char = ',';
const MISSING: char = '.';
const NUL: u8 = 0x00;

pub(super) fn write_values<W>(
    writer: &mut W,
    format: &Map<Format>,
    values: &[Option<Value<'_>>],
) -> io::Result<()>
where
    W: Write,
{
    use noodles_vcf::header::record::value::map::format::{self, Number};

    match format.ty() {
        format::Type::Integer => match format.number() {
            Number::Count(1) => write_integer_values(writer, values),
            _ => write_integer_array_values(writer, values),
        },
        format::Type::Float => match format.number() {
            Number::Count(1) => write_float_values(writer, values),
            _ => write_float_array_values(writer, values),
        },
        format::Type::Character => match format.number() {
            Number::Count(1) => write_character_values(writer, values),
            _ => write_character_array_values(writer, values),
        },
        format::Type::String => match format.number() {
            Number::Count(1) => write_string_values(writer, values),
            _ => write_string_array_values(writer, values),
        },
    }
}

fn write_integer_values<W>(writer: &mut W, values: &[Option<Value<'_>>]) -> io::Result<()>
where
    W: Write,
{
    let (mut min, mut max) = (i32::MAX, i32::MIN);

    for value in values {
        let n = match value {
            Some(Value::Integer(n)) => *n,
            _ => 0,
        };

        min = cmp::min(min, n);
        max = cmp::max(max, n);
    }

    if min >= i32::from(Int8::MIN_VALUE) {
        if max <= i32::from(Int8::MAX_VALUE) {
            write_int8_values(writer, values)
        } else if max <= i32::from(Int16::MAX_VALUE) {
            write_int16_values(writer, values)
        } else {
            write_int32_values(writer, values)
        }
    } else if min >= i32::from(Int16::MIN_VALUE) {
        if max <= i32::from(Int16::MAX_VALUE) {
            write_int16_values(writer, values)
        } else {
            write_int32_values(writer, values)
        }
    } else if min >= Int32::MIN_VALUE {
        write_int32_values(writer, values)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid genotype field integer value: {min}"),
        ))
    }
}

fn write_int8_values<W>(writer: &mut W, values: &[Option<Value<'_>>]) -> io::Result<()>
where
    W: Write,
{
    write_type(writer, Some(Type::Int8(1)))?;

    for value in values {
        match value {
            Some(Value::Integer(n)) => {
                let m =
                    i8::try_from(*n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                writer.write_i8(m)?;
            }
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Integer, got {v:?}"),
                ));
            }
            None => writer.write_i8(i8::from(Int8::Missing))?,
        }
    }

    Ok(())
}

fn write_int16_values<W>(writer: &mut W, values: &[Option<Value<'_>>]) -> io::Result<()>
where
    W: Write,
{
    write_type(writer, Some(Type::Int16(1)))?;

    for value in values {
        match value {
            Some(Value::Integer(n)) => {
                let m = i16::try_from(*n)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                writer.write_i16::<LittleEndian>(m)?;
            }
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Integer, got {v:?}"),
                ));
            }
            None => writer.write_i16::<LittleEndian>(i16::from(Int16::Missing))?,
        }
    }

    Ok(())
}

fn write_int32_values<W>(writer: &mut W, values: &[Option<Value<'_>>]) -> io::Result<()>
where
    W: Write,
{
    write_type(writer, Some(Type::Int32(1)))?;

    for value in values {
        match value {
            Some(Value::Integer(n)) => write_i32_le(writer, *n)?,
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Integer, got {v:?}"),
                ));
            }
            None => write_i32_le(writer, i32::from(Int32::Missing))?,
        }
    }

    Ok(())
}

fn write_integer_array_values<W>(writer: &mut W, values: &[Option<Value<'_>>]) -> io::Result<()>
where
    W: Write,
{
    let mut max_len = usize::MIN;
    let (mut min, mut max) = (i32::MAX, i32::MIN);

    for value in values {
        match value {
            Some(Value::Array(Array::Integer(vs))) => {
                max_len = cmp::max(max_len, vs.len());

                for result in vs.iter() {
                    let v = result?;
                    let n = v.unwrap_or_default();
                    min = cmp::min(min, n);
                    max = cmp::max(max, n);
                }
            }
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Array(Array::Integer), got {v:?}"),
                ));
            }
            None => {}
        }
    }

    if min >= i32::from(Int8::MIN_VALUE) {
        if max <= i32::from(Int8::MAX_VALUE) {
            write_int8_array_values(writer, values, max_len)
        } else if max <= i32::from(Int16::MAX_VALUE) {
            write_int16_array_values(writer, values, max_len)
        } else {
            write_int32_array_values(writer, values, max_len)
        }
    } else if min >= i32::from(Int16::MIN_VALUE) {
        if max <= i32::from(Int16::MAX_VALUE) {
            write_int16_array_values(writer, values, max_len)
        } else {
            write_int32_array_values(writer, values, max_len)
        }
    } else if min >= Int32::MIN_VALUE {
        write_int32_array_values(writer, values, max_len)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid genotype field integer array value: {min}"),
        ))
    }
}

fn write_int8_array_values<W>(
    writer: &mut W,
    values: &[Option<Value<'_>>],
    max_len: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_type(writer, Some(Type::Int8(max_len)))?;

    for value in values {
        let len = match value {
            Some(Value::Array(Array::Integer(vs))) => {
                for result in vs.iter() {
                    let v = result?;

                    let n = match v {
                        Some(n) => i8::try_from(n)
                            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?,
                        None => i8::from(Int8::Missing),
                    };

                    writer.write_i8(n)?;
                }

                vs.len()
            }
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Array(Array::Integer), got {v:?}"),
                ));
            }
            None => {
                writer.write_i8(i8::from(Int8::Missing))?;
                1
            }
        };

        if len < max_len {
            for _ in 0..(max_len - len) {
                writer.write_i8(i8::from(Int8::EndOfVector))?;
            }
        }
    }

    Ok(())
}

fn write_int16_array_values<W>(
    writer: &mut W,
    values: &[Option<Value<'_>>],
    max_len: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_type(writer, Some(Type::Int16(max_len)))?;

    for value in values {
        let len = match value {
            Some(Value::Array(Array::Integer(vs))) => {
                for result in vs.iter() {
                    let v = result?;

                    let n = match v {
                        Some(n) => i16::try_from(n)
                            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?,
                        None => i16::from(Int16::Missing),
                    };

                    writer.write_i16::<LittleEndian>(n)?;
                }

                vs.len()
            }
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Array(Array::Integer), got {v:?}"),
                ));
            }
            None => {
                writer.write_i16::<LittleEndian>(i16::from(Int16::Missing))?;
                1
            }
        };

        if len < max_len {
            for _ in 0..(max_len - len) {
                writer.write_i16::<LittleEndian>(i16::from(Int16::EndOfVector))?;
            }
        }
    }

    Ok(())
}

fn write_int32_array_values<W>(
    writer: &mut W,
    values: &[Option<Value<'_>>],
    max_len: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_type(writer, Some(Type::Int32(max_len)))?;

    for value in values {
        let len = match value {
            Some(Value::Array(Array::Integer(vs))) => {
                for result in vs.iter() {
                    let v = result?;

                    let n = match v {
                        Some(n) => n,
                        None => i32::from(Int32::Missing),
                    };

                    write_i32_le(writer, n)?;
                }

                vs.len()
            }
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Array(Array::Integer), got {v:?}"),
                ));
            }
            None => {
                write_i32_le(writer, i32::from(Int32::Missing))?;
                1
            }
        };

        if len < max_len {
            for _ in 0..(max_len - len) {
                write_i32_le(writer, i32::from(Int32::EndOfVector))?;
            }
        }
    }

    Ok(())
}

fn write_float_values<W>(writer: &mut W, values: &[Option<Value<'_>>]) -> io::Result<()>
where
    W: Write,
{
    write_type(writer, Some(Type::Float(1)))?;

    for value in values {
        match value {
            Some(Value::Float(n)) => {
                writer.write_f32::<LittleEndian>(*n)?;
            }
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Float, got {v:?}"),
                ));
            }
            None => writer.write_f32::<LittleEndian>(f32::from(Float::Missing))?,
        }
    }

    Ok(())
}

fn write_float_array_values<W>(writer: &mut W, values: &[Option<Value<'_>>]) -> io::Result<()>
where
    W: Write,
{
    let max_len = values
        .iter()
        .flat_map(|value| match value {
            Some(Value::Array(Array::Float(vs))) => Some(vs.len()),
            _ => None,
        })
        .max()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing float array values"))?;

    write_type(writer, Some(Type::Float(max_len)))?;

    for value in values {
        let len = match value {
            Some(Value::Array(Array::Float(vs))) => {
                for result in vs.iter() {
                    let v = result?;
                    let raw_value = v.unwrap_or(f32::from(Float::Missing));
                    writer.write_f32::<LittleEndian>(raw_value)?;
                }

                vs.len()
            }
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Array(Array::Float), got {v:?}"),
                ));
            }
            None => {
                writer.write_f32::<LittleEndian>(f32::from(Float::Missing))?;
                1
            }
        };

        if len < max_len {
            for _ in 0..(max_len - len) {
                writer.write_f32::<LittleEndian>(f32::from(Float::EndOfVector))?;
            }
        }
    }

    Ok(())
}

fn write_character_values<W>(writer: &mut W, values: &[Option<Value<'_>>]) -> io::Result<()>
where
    W: Write,
{
    let mut raw_strings = Vec::with_capacity(values.len());

    for value in values {
        match value {
            Some(Value::Character(c)) => {
                let s = String::from(*c);
                raw_strings.push(Some(s));
            }
            None => raw_strings.push(None),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Character, got {v:?}"),
                ));
            }
        }
    }

    let string_values: Vec<_> = raw_strings
        .iter()
        .map(|s| s.as_ref().map(|t| Value::String(Cow::from(t))))
        .collect();

    write_string_values(writer, &string_values)
}

fn write_character_array_values<W>(writer: &mut W, values: &[Option<Value<'_>>]) -> io::Result<()>
where
    W: Write,
{
    let mut raw_strings = Vec::with_capacity(values.len());

    for value in values {
        match value {
            Some(Value::Array(Array::Character(cs))) => {
                let mut s = String::new();

                for (i, result) in cs.iter().enumerate() {
                    if i > 0 {
                        s.push(DELIMITER);
                    }

                    match result? {
                        Some(c) => s.push(c),
                        None => s.push(MISSING),
                    }
                }

                raw_strings.push(Some(s));
            }
            None => raw_strings.push(None),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Character, got {v:?}"),
                ));
            }
        }
    }

    let string_values: Vec<_> = raw_strings
        .iter()
        .map(|s| s.as_ref().map(|t| Value::String(Cow::from(t))))
        .collect();

    write_string_values(writer, &string_values)
}

fn write_string_values<W>(writer: &mut W, values: &[Option<Value<'_>>]) -> io::Result<()>
where
    W: Write,
{
    let max_len = values
        .iter()
        .flat_map(|value| match value {
            Some(Value::String(s)) => Some(s.len()),
            _ => None,
        })
        .max()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing String values"))?;

    let mut buf = Vec::with_capacity(values.len() * max_len);

    for value in values {
        match value {
            Some(Value::String(s)) => {
                buf.extend(s.bytes());

                if s.len() < max_len {
                    buf.resize(buf.len() + (max_len - s.len()), NUL);
                }
            }
            None => {
                buf.push(b'.');

                if 1 < max_len {
                    buf.resize(buf.len() + (max_len - 1), NUL);
                }
            }
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected String, got {v:?}"),
                ));
            }
        }
    }

    write_type(writer, Some(Type::String(max_len)))?;
    writer.write_all(&buf)?;

    Ok(())
}

fn write_string_array_values<W>(writer: &mut W, values: &[Option<Value<'_>>]) -> io::Result<()>
where
    W: Write,
{
    let mut serialized_values = Vec::with_capacity(values.len());

    for value in values {
        match value {
            Some(Value::Array(Array::String(vs))) => {
                let mut s = String::new();

                for (i, result) in vs.iter().enumerate() {
                    if i > 0 {
                        s.push(DELIMITER);
                    }

                    match result? {
                        Some(t) => s.push_str(&t),
                        None => s.push(MISSING),
                    }
                }

                serialized_values.push(s);
            }
            None => serialized_values.push(String::from(MISSING)),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected String, got {v:?}"),
                ));
            }
        }
    }

    let max_len = serialized_values
        .iter()
        .map(|s| s.len())
        .max()
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "missing serialized string array values",
            )
        })?;

    let mut buf = Vec::with_capacity(serialized_values.len() * max_len);

    for s in serialized_values {
        buf.extend(s.bytes());

        if s.len() < max_len {
            buf.resize(buf.len() + (max_len - s.len()), NUL);
        }
    }

    write_type(writer, Some(Type::String(max_len)))?;
    writer.write_all(&buf)?;

    Ok(())
}

pub(super) fn write_genotype_values<W>(
    writer: &mut W,
    values: &[Option<Value<'_>>],
) -> io::Result<()>
where
    W: Write,
{
    let mut raw_values = Vec::with_capacity(values.len());
    let mut max_len = 0;

    for value in values {
        let raw_value = match value {
            Some(Value::String(s)) => encode_genotype_str(s)?,
            Some(Value::Genotype(genotype)) => encode_genotype(genotype.as_ref())?,
            _ => return Err(io::Error::from(io::ErrorKind::InvalidInput)),
        };

        max_len = cmp::max(max_len, raw_value.len());

        raw_values.push(raw_value);
    }

    write_type(writer, Some(Type::Int8(max_len)))?;

    for raw_value in raw_values {
        let len = raw_value.len();
        let pad = max_len - len;

        for n in raw_value {
            let m = u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
            writer.write_all(&[m])?;

            for _ in 0..pad {
                writer.write_all(&[i8::from(Int8::EndOfVector) as u8])?;
            }
        }
    }

    Ok(())
}

fn encode_genotype_str(genotype: &str) -> io::Result<Vec<i8>> {
    const MISSING_ALLELE: &str = ".";

    fn is_phasing(c: char) -> bool {
        matches!(c, '|' | '/')
    }

    fn encode(s: &str, phasing: &str) -> io::Result<i8> {
        if s == MISSING_ALLELE {
            return Ok(0);
        }

        let j: i8 = s
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let is_phased = phasing == "|";

        let mut i = (j + 1) << 1;

        if is_phased {
            i |= 0x01;
        }

        Ok(i)
    }

    let mut values = Vec::new();
    let mut start = 0;
    let mut last_phasing = "/";

    for (end, phasing) in genotype.match_indices(is_phasing) {
        let i = encode(&genotype[start..end], last_phasing)?;
        values.push(i);
        start = end + 1;
        last_phasing = phasing;
    }

    let i = encode(&genotype[start..], last_phasing)?;
    values.push(i);

    Ok(values)
}

fn encode_genotype(genotype: &dyn Genotype) -> io::Result<Vec<i8>> {
    fn encode(position: Option<usize>, phasing: Phasing) -> io::Result<i8> {
        let i = if let Some(position) = position {
            i8::try_from(position).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
        } else {
            return Ok(0);
        };

        let mut n = (i + 1) << 1;

        if phasing == Phasing::Phased {
            n |= 0x01;
        }

        Ok(n)
    }

    genotype
        .iter()
        .map(|result| result.and_then(|(position, phasing)| encode(position, phasing)))
        .collect()
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{
        header::record::value::map::format::{self, Number},
        variant::record_buf::samples::sample::Value as ValueBuf,
    };

    use super::*;

    #[test]
    fn test_write_values_with_integer_values() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            buf: &mut Vec<u8>,
            format: &Map<Format>,
            values: &[Option<Value<'_>>],
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_values(buf, format, values)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let key = Map::<Format>::new(Number::Count(1), format::Type::Integer, String::new());

        let mut buf = Vec::new();

        let values = [
            Some(Value::Integer(-2147483641)),
            Some(Value::Integer(-2147483640)),
            None,
        ];
        buf.clear();
        assert!(matches!(
            write_values(&mut buf, &key, &values),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        let values = [
            Some(Value::Integer(-2147483640)),
            Some(Value::Integer(-2147483639)),
            None,
        ];
        let expected = [
            0x13, // Some(Type::Int32(1))
            0x08, 0x00, 0x00, 0x80, // Some(-2147483640)
            0x09, 0x00, 0x00, 0x80, // Some(-2147483639)
            0x00, 0x00, 0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [
            Some(Value::Integer(-32761)),
            Some(Value::Integer(-32760)),
            None,
        ];
        let expected = [
            0x13, // Some(Type::Int32(1))
            0x07, 0x80, 0xff, 0xff, // Some(-32761)
            0x08, 0x80, 0xff, 0xff, // Some(-32760)
            0x00, 0x00, 0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [
            Some(Value::Integer(-32760)),
            Some(Value::Integer(-32759)),
            None,
        ];
        let expected = [
            0x12, // Some(Type::Int16(1))
            0x08, 0x80, // Some(-32760)
            0x09, 0x80, // Some(-32759)
            0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [Some(Value::Integer(-121)), Some(Value::Integer(-120)), None];
        let expected = [
            0x12, // Some(Type::Int16(1))
            0x87, 0xff, // Some(-121)
            0x88, 0xff, // Some(-120)
            0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [Some(Value::Integer(-120)), Some(Value::Integer(-119)), None];
        let expected = [
            0x11, // Some(Type::Int8(1))
            0x88, // Some(-120)
            0x89, // Some(-119)
            0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [
            Some(Value::Integer(-1)),
            Some(Value::Integer(0)),
            Some(Value::Integer(1)),
            None,
        ];
        let expected = [
            0x11, // Some(Type::Int8(1))
            0xff, // Some(-1)
            0x00, // Some(0)
            0x01, // Some(1)
            0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [Some(Value::Integer(126)), Some(Value::Integer(127)), None];
        let expected = [
            0x11, // Some(Type::Int8(1))
            0x7e, // Some(126)
            0x7f, // Some(127)
            0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [Some(Value::Integer(127)), Some(Value::Integer(128)), None];
        let expected = [
            0x12, // Some(Type::Int16(1))
            0x7f, 0x00, // Some(127)
            0x80, 0x00, // Some(128)
            0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [
            Some(Value::Integer(32766)),
            Some(Value::Integer(32767)),
            None,
        ];
        let expected = [
            0x12, // Some(Type::Int16(1))
            0xfe, 0x7f, // Some(32766)
            0xff, 0x7f, // Some(32767)
            0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [
            Some(Value::Integer(32767)),
            Some(Value::Integer(32768)),
            None,
        ];
        let expected = [
            0x13, // Some(Type::Int32(1))
            0xff, 0x7f, 0x00, 0x00, // Some(32767)
            0x00, 0x80, 0x00, 0x00, // Some(32768)
            0x00, 0x00, 0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [
            Some(Value::Integer(2147483646)),
            Some(Value::Integer(2147483647)),
            None,
        ];
        let expected = [
            0x13, // Some(Type::Int32(1))
            0xfe, 0xff, 0xff, 0x7f, // Some(2147483646)
            0xff, 0xff, 0xff, 0x7f, // Some(2147483647)
            0x00, 0x00, 0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        Ok(())
    }

    #[test]
    fn test_write_values_with_integer_array_values() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            buf: &mut Vec<u8>,
            format: &Map<Format>,
            values: &[Option<Value>],
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_values(buf, format, values)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let format = Map::<Format>::new(Number::Count(2), format::Type::Integer, String::new());

        let mut buf = Vec::new();

        let value_0 = ValueBuf::from(vec![Some(-2147483641), Some(-2147483640)]);
        let value_1 = ValueBuf::from(vec![Some(-2147483639), None]);
        let value_2 = ValueBuf::from(vec![Some(-2147483638)]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];
        buf.clear();
        assert!(matches!(
            write_values(&mut buf, &format, &values),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        let value_0 = ValueBuf::from(vec![Some(-2147483640), Some(-2147483639)]);
        let value_1 = ValueBuf::from(vec![Some(-2147483638), None]);
        let value_2 = ValueBuf::from(vec![Some(-2147483637)]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];
        #[rustfmt::skip]
        let expected = [
            0x23, // Some(Type::Int32(2))
            0x08, 0x00, 0x00, 0x80, 0x09, 0x00, 0x00, 0x80, // Some([Some(-2147483640), Some(-2147483639)])
            0x0a, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, // Some([Some(-2147483638), None])
            0x0b, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // Some([Some(-2147483637)])
            0x00, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // None
        ];
        t(&mut buf, &format, &values, &expected)?;

        let value_0 = ValueBuf::from(vec![Some(-32761), Some(-32760)]);
        let value_1 = ValueBuf::from(vec![Some(-32759), None]);
        let value_2 = ValueBuf::from(vec![Some(-32758)]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];
        #[rustfmt::skip]
        let expected = [
            0x23, // Some(Type::Int32(2))
            0x07, 0x80, 0xff, 0xff, 0x08, 0x80, 0xff, 0xff, // Some([Some(-32761), Some(-32760)])
            0x09, 0x80, 0xff, 0xff, 0x00, 0x00, 0x00, 0x80, // Some([Some(-32759), None])
            0x0a, 0x80, 0xff, 0xff, 0x01, 0x00, 0x00, 0x80, // Some([Some(-32761)])
            0x00, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // None
        ];
        t(&mut buf, &format, &values, &expected)?;

        let value_0 = ValueBuf::from(vec![Some(-32760), Some(-32759)]);
        let value_1 = ValueBuf::from(vec![Some(-32758), None]);
        let value_2 = ValueBuf::from(vec![Some(-32757)]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];
        let expected = [
            0x22, // Some(Type::Int16(2))
            0x08, 0x80, 0x09, 0x80, // Some([Some(-32760), Some(-32759)])
            0x0a, 0x80, 0x00, 0x80, // Some([Some(-32758), None])
            0x0b, 0x80, 0x01, 0x80, // Some([Some(-32757)])
            0x00, 0x80, 0x01, 0x80, // None
        ];
        t(&mut buf, &format, &values, &expected)?;

        let value_0 = ValueBuf::from(vec![Some(-121), Some(-120)]);
        let value_1 = ValueBuf::from(vec![Some(-119), None]);
        let value_2 = ValueBuf::from(vec![Some(-118)]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];
        let expected = [
            0x22, // Some(Type::Int16(2))
            0x87, 0xff, 0x88, 0xff, // Some([Some(-121), Some(-120)])
            0x89, 0xff, 0x00, 0x80, // Some([Some(-119), None])
            0x8a, 0xff, 0x01, 0x80, // Some([Some(-118)])
            0x00, 0x80, 0x01, 0x80, // None
        ];
        t(&mut buf, &format, &values, &expected)?;

        let value_0 = ValueBuf::from(vec![Some(-120), Some(-119)]);
        let value_1 = ValueBuf::from(vec![Some(-118), None]);
        let value_2 = ValueBuf::from(vec![Some(-117)]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];
        let expected = [
            0x21, // Some(Type::Int8(2))
            0x88, 0x89, // Some([Some(-120), Some(-119)])
            0x8a, 0x80, // Some([Some(-118), None])
            0x8b, 0x81, // Some([Some(-117)])
            0x80, 0x81, // None
        ];
        t(&mut buf, &format, &values, &expected)?;

        let value_0 = ValueBuf::from(vec![Some(-1), Some(0)]);
        let value_1 = ValueBuf::from(vec![Some(1), None]);
        let value_2 = ValueBuf::from(vec![Some(2)]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];
        let expected = [
            0x21, // Some(Type::Int8(2))
            0xff, 0x00, // Some([Some(-1), Some(0)])
            0x01, 0x80, // Some([Some(1), None])
            0x02, 0x81, // Some([Some(2)])
            0x80, 0x81, // None
        ];
        t(&mut buf, &format, &values, &expected)?;

        let value_0 = ValueBuf::from(vec![Some(124), Some(125)]);
        let value_1 = ValueBuf::from(vec![Some(126), None]);
        let value_2 = ValueBuf::from(vec![Some(127)]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];
        let expected = [
            0x21, // Some(Type::Int8(2))
            0x7c, 0x7d, // Some([Some(124), Some(125)])
            0x7e, 0x80, // Some([Some(126), None])
            0x7f, 0x81, // Some([Some(127)])
            0x80, 0x81, // None
        ];
        t(&mut buf, &format, &values, &expected)?;

        let value_0 = ValueBuf::from(vec![Some(125), Some(126)]);
        let value_1 = ValueBuf::from(vec![Some(127), None]);
        let value_2 = ValueBuf::from(vec![Some(128)]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];
        let expected = [
            0x22, // Some(Type::Int16(2))
            0x7d, 0x00, 0x7e, 0x00, // Some([Some(125), Some(126)])
            0x7f, 0x00, 0x00, 0x80, // Some([Some(127), None])
            0x80, 0x00, 0x01, 0x80, // Some([Some(128)])
            0x00, 0x80, 0x01, 0x80, // None
        ];
        t(&mut buf, &format, &values, &expected)?;

        let value_0 = ValueBuf::from(vec![Some(32764), Some(32765)]);
        let value_1 = ValueBuf::from(vec![Some(32766), None]);
        let value_2 = ValueBuf::from(vec![Some(32767)]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];
        let expected = [
            0x22, // Some(Type::Int16(2))
            0xfc, 0x7f, 0xfd, 0x7f, // Some([Some(32764), Some(32765)])
            0xfe, 0x7f, 0x00, 0x80, // Some([Some(32766), None])
            0xff, 0x7f, 0x01, 0x80, // Some([Some(32767)])
            0x00, 0x80, 0x01, 0x80, // None
        ];
        t(&mut buf, &format, &values, &expected)?;

        let value_0 = ValueBuf::from(vec![Some(32765), Some(32766)]);
        let value_1 = ValueBuf::from(vec![Some(32767), None]);
        let value_2 = ValueBuf::from(vec![Some(32768)]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];
        #[rustfmt::skip]
        let expected = [
            0x23, // Some(Type::Int32(2))
            0xfd, 0x7f, 0x00, 0x00, 0xfe, 0x7f, 0x00, 0x00, // Some([Some(32765), Some(32766)])
            0xff, 0x7f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, // Some([Some(32767), None])
            0x00, 0x80, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, // Some([Some(32768)])
            0x00, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // None
        ];
        t(&mut buf, &format, &values, &expected)?;

        let value_0 = ValueBuf::from(vec![Some(2147483644), Some(2147483645)]);
        let value_1 = ValueBuf::from(vec![Some(2147483646), None]);
        let value_2 = ValueBuf::from(vec![Some(2147483647)]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];
        #[rustfmt::skip]
        let expected = [
            0x23, // Some(Type::Int32(2))
            0xfc, 0xff, 0xff, 0x7f, 0xfd, 0xff, 0xff, 0x7f, // Some([Some(2147483645), Some(2147483645)])
            0xfe, 0xff, 0xff, 0x7f, 0x00, 0x00, 0x00, 0x80, // Some([Some(2147483646), None])
            0xff, 0xff, 0xff, 0x7f, 0x01, 0x00, 0x00, 0x80, // Some([Some(2147483647)])
            0x00, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // None
        ];
        t(&mut buf, &format, &values, &expected)?;

        Ok(())
    }

    #[test]
    fn test_write_values_with_float_values() -> Result<(), Box<dyn std::error::Error>> {
        let format = Map::<Format>::new(Number::Count(1), format::Type::Float, String::new());

        let values = [Some(Value::Float(0.0)), Some(Value::Float(1.0)), None];

        let mut buf = Vec::new();
        write_values(&mut buf, &format, &values)?;

        let expected = [
            0x15, // Some(Type::Float(1))
            0x00, 0x00, 0x00, 0x00, // Some(0.0)
            0x00, 0x00, 0x80, 0x3f, // Some(1.0)
            0x01, 0x00, 0x80, 0x7f, // None
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_values_with_float_array_values() -> Result<(), Box<dyn std::error::Error>> {
        let format = Map::<Format>::new(Number::Count(2), format::Type::Float, String::new());

        let value_0 = ValueBuf::from(vec![Some(0.0), Some(1.0)]);
        let value_1 = ValueBuf::from(vec![Some(0.0), None]);
        let value_2 = ValueBuf::from(vec![Some(0.0)]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];

        let mut buf = Vec::new();
        write_values(&mut buf, &format, &values)?;

        let expected = [
            0x25, // Some(Type::Float(2))
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x3f, // Some([Some(0.0), Some(1.0)])
            0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x80, 0x7f, // Some([Some(0.0), None])
            0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x80, 0x7f, // Some([Some(0.0)])
            0x01, 0x00, 0x80, 0x7f, 0x02, 0x00, 0x80, 0x7f, // None
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_values_with_character_values() -> Result<(), Box<dyn std::error::Error>> {
        let format = Map::<Format>::new(Number::Count(1), format::Type::Character, String::new());

        let values = [
            Some(Value::Character('n')),
            Some(Value::Character('d')),
            None,
        ];

        let mut buf = Vec::new();
        write_values(&mut buf, &format, &values)?;

        let expected = [
            0x17, // Some(Type::String(1))
            0x6e, // Some('n')
            0x64, // Some('d')
            0x2e, // None
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_values_with_character_array_values() -> Result<(), Box<dyn std::error::Error>> {
        let format = Map::<Format>::new(Number::Count(2), format::Type::Character, String::new());

        let value_0 = ValueBuf::from(vec![Some('n'), Some('d')]);
        let value_1 = ValueBuf::from(vec![Some('l'), None]);
        let value_2 = ValueBuf::from(vec![Some('s')]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];

        let mut buf = Vec::new();
        write_values(&mut buf, &format, &values)?;

        let expected = [
            0x37, // Some(Type::String(1))
            b'n', b',', b'd', // Some([Some('n'), Some('d')])
            b'l', b',', b'.', // Some([Some('d'), None])
            b's', 0x00, 0x00, // Some([Some('s')])
            b'.', 0x00, 0x00, // None
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_values_with_string_values() -> Result<(), Box<dyn std::error::Error>> {
        let format = Map::<Format>::new(Number::Count(1), format::Type::String, String::new());

        let value_0 = ValueBuf::from("n");
        let value_1 = ValueBuf::from("ndl");
        let value_2 = ValueBuf::from("ndls");
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];

        let mut buf = Vec::new();
        write_values(&mut buf, &format, &values)?;

        let expected = [
            0x47, // Some(Type::String(4))
            b'n', 0x00, 0x00, 0x00, // Some("n")
            b'n', b'd', b'l', 0x00, // Some("ndl")
            b'n', b'd', b'l', b's', // Some("ndls")
            b'.', 0x00, 0x00, 0x00, // None
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_values_with_string_array_values() -> Result<(), Box<dyn std::error::Error>> {
        let format = Map::<Format>::new(Number::Count(2), format::Type::String, String::new());

        let value_0 = ValueBuf::from(vec![Some(String::from("n")), Some(String::from("nd"))]);
        let value_1 = ValueBuf::from(vec![Some(String::from("ndls")), None]);
        let value_2 = ValueBuf::from(vec![Some(String::from("nd"))]);
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
            None,
        ];

        let mut buf = Vec::new();
        write_values(&mut buf, &format, &values)?;

        let expected = [
            0x67, // Some(Type::String(6))
            b'n', b',', b'n', b'd', 0x00, 0x00, // Some([Some("n"), Some("nd")])
            b'n', b'd', b'l', b's', b',', b'.', // Some([Some("ndls"), None])
            b'n', b'd', 0x00, 0x00, 0x00, 0x00, // Some([Some("nd")])
            b'.', 0x00, 0x00, 0x00, 0x00, 0x00, // None
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_genotype_values() -> io::Result<()> {
        let value_0 = ValueBuf::from("0/1");
        let value_1 = ValueBuf::from("1/1");
        let value_2 = ValueBuf::from("0/0");
        let values = [
            Some((&value_0).into()),
            Some((&value_1).into()),
            Some((&value_2).into()),
        ];

        let mut buf = Vec::new();
        write_genotype_values(&mut buf, &values)?;

        let expected = [
            0x21, // Some(Type::Int8(2))
            0x02, 0x04, // "0/1"
            0x04, 0x04, // "1/1"
            0x02, 0x02, // "0/0"
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_genotype_values_with_padding() -> io::Result<()> {
        let value_0 = ValueBuf::from("0");
        let value_1 = ValueBuf::from("0/1");
        let values = [Some((&value_0).into()), Some((&value_1).into())];

        let mut buf = Vec::new();
        write_genotype_values(&mut buf, &values)?;

        let expected = [
            0x21, // Some(Type::Int8(2))
            0x02, 0x81, // "0"
            0x02, 0x04, // "0/1"
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_encode_genotype_str() -> io::Result<()> {
        assert_eq!(encode_genotype_str("0/1")?, [0x02, 0x04]);
        assert_eq!(encode_genotype_str("0|1")?, [0x02, 0x05]);
        assert_eq!(encode_genotype_str("./.")?, [0x00, 0x00]);
        assert_eq!(encode_genotype_str("0")?, [0x02]);
        assert_eq!(encode_genotype_str("1")?, [0x04]);
        assert_eq!(encode_genotype_str("0/1/2")?, [0x02, 0x04, 0x06]);
        assert_eq!(encode_genotype_str("0/1|2")?, [0x02, 0x04, 0x07]);

        Ok(())
    }

    #[test]
    fn test_encode_genotype() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_vcf::variant::{
            record::samples::series::value::genotype::Phasing,
            record_buf::samples::sample::value::genotype::Allele,
        };

        let genotype = &[
            Allele::new(Some(0), Phasing::Unphased),
            Allele::new(Some(1), Phasing::Unphased),
        ]
        .into_iter()
        .collect();
        assert_eq!(encode_genotype(&genotype)?, [0x02, 0x04]);

        let genotype = &[
            Allele::new(Some(0), Phasing::Phased),
            Allele::new(Some(1), Phasing::Phased),
        ]
        .into_iter()
        .collect();
        assert_eq!(encode_genotype(&genotype)?, [0x03, 0x05]);

        let genotype = &[
            Allele::new(None, Phasing::Unphased),
            Allele::new(None, Phasing::Unphased),
        ]
        .into_iter()
        .collect();
        assert_eq!(encode_genotype(&genotype)?, [0x00, 0x00]);

        let genotype = &[Allele::new(Some(0), Phasing::Phased)]
            .into_iter()
            .collect();
        assert_eq!(encode_genotype(&genotype)?, [0x03]);

        let genotype = &[Allele::new(Some(1), Phasing::Phased)]
            .into_iter()
            .collect();
        assert_eq!(encode_genotype(&genotype)?, [0x05]);

        let genotype = &[
            Allele::new(Some(0), Phasing::Unphased),
            Allele::new(Some(1), Phasing::Unphased),
            Allele::new(Some(2), Phasing::Unphased),
        ]
        .into_iter()
        .collect();
        assert_eq!(encode_genotype(&genotype)?, [0x02, 0x04, 0x06]);

        let genotype = &[
            Allele::new(Some(0), Phasing::Unphased),
            Allele::new(Some(1), Phasing::Unphased),
            Allele::new(Some(2), Phasing::Phased),
        ]
        .into_iter()
        .collect();
        assert_eq!(encode_genotype(&genotype)?, [0x02, 0x04, 0x07]);

        Ok(())
    }
}
