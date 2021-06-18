use std::{
    cmp,
    convert::TryFrom,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_vcf::{
    self as vcf,
    record::genotype::field::{Key, Value},
};

use crate::{
    header::StringMap,
    record::value::{Float, Int16, Int32, Int8, Type},
    writer::{string_map::write_string_map_index, value::write_type},
};

const NUL: u8 = 0x00;
const DELIMITER: char = ',';
const MISSING_VALUE: char = '.';

pub fn write_genotypes<W>(
    writer: &mut W,
    string_map: &StringMap,
    formats: &vcf::record::Format,
    genotypes: &[vcf::record::Genotype],
) -> io::Result<()>
where
    W: Write,
{
    for (i, key) in formats.iter().enumerate() {
        write_genotype_field_key(writer, string_map, key)?;

        let mut values = Vec::with_capacity(formats.len());

        for genotype in genotypes {
            let value = genotype
                .get_index(i)
                .map(|(_, field)| field.value())
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "missing genotype field")
                })?;

            values.push(value);
        }

        write_genotype_field_values(writer, key, &values)?;
    }

    Ok(())
}

pub fn write_genotype_field_key<W>(
    writer: &mut W,
    string_map: &StringMap,
    key: &Key,
) -> io::Result<()>
where
    W: Write,
{
    string_map
        .get_index_of(key.as_ref())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("genotype key not in string map: {:?}", key),
            )
        })
        .and_then(|i| write_string_map_index(writer, i))
}

pub fn write_genotype_field_values<W>(
    writer: &mut W,
    key: &Key,
    values: &[Option<&Value>],
) -> io::Result<()>
where
    W: Write,
{
    use vcf::header::{format, Number};

    match key.ty() {
        format::Type::Integer => match key.number() {
            Number::Count(1) => write_genotype_field_integer_values(writer, values),
            _ => write_genotype_field_integer_array_values(writer, values),
        },
        format::Type::Float => match key.number() {
            Number::Count(1) => write_genotype_field_float_values(writer, values),
            _ => write_genotype_field_float_array_values(writer, values),
        },
        format::Type::Character => match key.number() {
            Number::Count(1) => write_genotype_field_character_values(writer, values),
            _ => write_genotype_field_character_array_values(writer, values),
        },
        format::Type::String => match key.number() {
            Number::Count(1) => write_genotype_field_string_values(writer, values),
            _ => write_genotype_field_string_array_values(writer, values),
        },
    }
}

fn write_genotype_field_integer_values<W>(
    writer: &mut W,
    values: &[Option<&Value>],
) -> io::Result<()>
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
            write_genotype_field_int8_values(writer, values)
        } else if max <= i32::from(Int16::MAX_VALUE) {
            write_genotype_field_int16_values(writer, values)
        } else {
            write_genotype_field_int32_values(writer, values)
        }
    } else if min >= i32::from(Int16::MIN_VALUE) {
        if max <= i32::from(Int16::MAX_VALUE) {
            write_genotype_field_int16_values(writer, values)
        } else {
            write_genotype_field_int32_values(writer, values)
        }
    } else if min >= Int32::MIN_VALUE {
        write_genotype_field_int32_values(writer, values)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid genotype field integer value: {}", min),
        ))
    }
}

fn write_genotype_field_int8_values<W>(writer: &mut W, values: &[Option<&Value>]) -> io::Result<()>
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
                    format!("type mismatch: expected Integer, got {:?}", v),
                ));
            }
            None => writer.write_i8(i8::from(Int8::Missing))?,
        }
    }

    Ok(())
}

fn write_genotype_field_int16_values<W>(writer: &mut W, values: &[Option<&Value>]) -> io::Result<()>
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
                    format!("type mismatch: expected Integer, got {:?}", v),
                ));
            }
            None => writer.write_i16::<LittleEndian>(i16::from(Int16::Missing))?,
        }
    }

    Ok(())
}

fn write_genotype_field_int32_values<W>(writer: &mut W, values: &[Option<&Value>]) -> io::Result<()>
where
    W: Write,
{
    write_type(writer, Some(Type::Int32(1)))?;

    for value in values {
        match value {
            Some(Value::Integer(n)) => {
                writer.write_i32::<LittleEndian>(*n)?;
            }
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Integer, got {:?}", v),
                ));
            }
            None => writer.write_i32::<LittleEndian>(i32::from(Int32::Missing))?,
        }
    }

    Ok(())
}

fn write_genotype_field_integer_array_values<W>(
    writer: &mut W,
    values: &[Option<&Value>],
) -> io::Result<()>
where
    W: Write,
{
    let mut max_len = usize::MIN;
    let (mut min, mut max) = (i32::MAX, i32::MIN);

    for value in values {
        match value {
            Some(Value::IntegerArray(vs)) => {
                max_len = cmp::max(max_len, vs.len());

                for v in vs {
                    let n = v.unwrap_or_default();
                    min = cmp::min(min, n);
                    max = cmp::max(max, n);
                }
            }
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected IntegerArray, got {:?}", v),
                ));
            }
            None => {}
        }
    }

    if min >= i32::from(Int8::MIN_VALUE) {
        if max <= i32::from(Int8::MAX_VALUE) {
            write_genotype_field_int8_array_values(writer, values, max_len)
        } else if max <= i32::from(Int16::MAX_VALUE) {
            write_genotype_field_int16_array_values(writer, values, max_len)
        } else {
            write_genotype_field_int32_array_values(writer, values, max_len)
        }
    } else if min >= i32::from(Int16::MIN_VALUE) {
        if max <= i32::from(Int16::MAX_VALUE) {
            write_genotype_field_int16_array_values(writer, values, max_len)
        } else {
            write_genotype_field_int32_array_values(writer, values, max_len)
        }
    } else if min >= Int32::MIN_VALUE {
        write_genotype_field_int32_array_values(writer, values, max_len)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid genotype field integer array value: {}", min),
        ))
    }
}

fn write_genotype_field_int8_array_values<W>(
    writer: &mut W,
    values: &[Option<&Value>],
    max_len: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_type(writer, Some(Type::Int8(max_len)))?;

    for value in values {
        let len = match value {
            Some(Value::IntegerArray(vs)) => {
                for v in vs {
                    let n = match v {
                        Some(n) => i8::try_from(*n)
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
                    format!("type mismatch: expected IntegerArray, got {:?}", v),
                ))
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

fn write_genotype_field_int16_array_values<W>(
    writer: &mut W,
    values: &[Option<&Value>],
    max_len: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_type(writer, Some(Type::Int16(max_len)))?;

    for value in values {
        let len = match value {
            Some(Value::IntegerArray(vs)) => {
                for v in vs {
                    let n = match v {
                        Some(n) => i16::try_from(*n)
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
                    format!("type mismatch: expected IntegerArray, got {:?}", v),
                ))
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

fn write_genotype_field_int32_array_values<W>(
    writer: &mut W,
    values: &[Option<&Value>],
    max_len: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_type(writer, Some(Type::Int32(max_len)))?;

    for value in values {
        let len = match value {
            Some(Value::IntegerArray(vs)) => {
                for v in vs {
                    let n = match v {
                        Some(n) => *n,
                        None => i32::from(Int32::Missing),
                    };

                    writer.write_i32::<LittleEndian>(n)?;
                }

                vs.len()
            }
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected IntegerArray, got {:?}", v),
                ))
            }
            None => {
                writer.write_i32::<LittleEndian>(i32::from(Int32::Missing))?;
                1
            }
        };

        if len < max_len {
            for _ in 0..(max_len - len) {
                writer.write_i32::<LittleEndian>(i32::from(Int32::EndOfVector))?;
            }
        }
    }

    Ok(())
}

fn write_genotype_field_float_values<W>(writer: &mut W, values: &[Option<&Value>]) -> io::Result<()>
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
                    format!("type mismatch: expected Float, got {:?}", v),
                ))
            }
            None => writer.write_f32::<LittleEndian>(f32::from(Float::Missing))?,
        }
    }

    Ok(())
}

fn write_genotype_field_float_array_values<W>(
    writer: &mut W,
    values: &[Option<&Value>],
) -> io::Result<()>
where
    W: Write,
{
    let max_len = values
        .iter()
        .flat_map(|value| match value {
            Some(Value::FloatArray(vs)) => Some(vs.len()),
            _ => None,
        })
        .max()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing FloatArray values"))?;

    write_type(writer, Some(Type::Float(max_len)))?;

    for value in values {
        let len = match value {
            Some(Value::FloatArray(vs)) => {
                for v in vs {
                    let raw_value = v.unwrap_or(f32::from(Float::Missing));
                    writer.write_f32::<LittleEndian>(raw_value)?;
                }

                vs.len()
            }
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected FloatArray, got {:?}", v),
                ))
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

fn write_genotype_field_character_values<W>(
    writer: &mut W,
    values: &[Option<&Value>],
) -> io::Result<()>
where
    W: Write,
{
    let mut string_values = Vec::with_capacity(values.len());

    for value in values {
        match value {
            Some(Value::Character(c)) => {
                let string_value = Value::String(String::from(*c));
                string_values.push(Some(string_value));
            }
            None => string_values.push(None),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Character, got {:?}", v),
                ))
            }
        }
    }

    let string_values_as_ref: Vec<_> = string_values.iter().map(|v| v.as_ref()).collect();

    write_genotype_field_string_values(writer, &string_values_as_ref)
}

fn write_genotype_field_character_array_values<W>(
    writer: &mut W,
    values: &[Option<&Value>],
) -> io::Result<()>
where
    W: Write,
{
    let mut string_values = Vec::with_capacity(values.len());

    for value in values {
        match value {
            Some(Value::CharacterArray(cs)) => {
                let mut s = String::new();

                for (i, c) in cs.iter().enumerate() {
                    if i > 0 {
                        s.push(DELIMITER);
                    }

                    match c {
                        Some(d) => s.push(*d),
                        None => s.push(MISSING_VALUE),
                    }
                }

                let string_value = Value::String(s);
                string_values.push(Some(string_value));
            }
            None => string_values.push(None),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Character, got {:?}", v),
                ))
            }
        }
    }

    let string_values_as_ref: Vec<_> = string_values.iter().map(|v| v.as_ref()).collect();

    write_genotype_field_string_values(writer, &string_values_as_ref)
}

fn write_genotype_field_string_values<W>(
    writer: &mut W,
    values: &[Option<&Value>],
) -> io::Result<()>
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
                    format!("type mismatch: expected String, got {:?}", v),
                ))
            }
        }
    }

    write_type(writer, Some(Type::String(max_len)))?;
    writer.write_all(&buf)?;

    Ok(())
}

fn write_genotype_field_string_array_values<W>(
    writer: &mut W,
    values: &[Option<&Value>],
) -> io::Result<()>
where
    W: Write,
{
    let mut serialized_values = Vec::with_capacity(values.len());

    for value in values {
        match value {
            Some(Value::StringArray(vs)) => {
                let mut s = String::new();

                for (i, v) in vs.iter().enumerate() {
                    if i > 0 {
                        s.push(DELIMITER);
                    }

                    match v {
                        Some(t) => s.push_str(t),
                        None => s.push(MISSING_VALUE),
                    }
                }

                serialized_values.push(s);
            }
            None => serialized_values.push(String::from(MISSING_VALUE)),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected String, got {:?}", v),
                ))
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
                "missing serialized StringArray values",
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

#[cfg(test)]
mod tests {
    use noodles_vcf::header::{format, Number};

    use super::*;

    #[test]
    fn test_write_genotype_field_values_with_integer_values() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            key: &Key,
            values: &[Option<&Value>],
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_genotype_field_values(buf, &key, &values)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let key = Key::Other(
            String::from("I32"),
            Number::Count(1),
            format::Type::Integer,
            String::default(),
        );

        let mut buf = Vec::new();

        let values = [
            Some(&Value::Integer(-2147483641)),
            Some(&Value::Integer(-2147483640)),
            None,
        ];
        buf.clear();
        assert!(matches!(
            write_genotype_field_values(&mut buf, &key, &values),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        let values = [
            Some(&Value::Integer(-2147483640)),
            Some(&Value::Integer(-2147483639)),
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
            Some(&Value::Integer(-32761)),
            Some(&Value::Integer(-32760)),
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
            Some(&Value::Integer(-32760)),
            Some(&Value::Integer(-32759)),
            None,
        ];
        let expected = [
            0x12, // Some(Type::Int16(1))
            0x08, 0x80, // Some(-32760)
            0x09, 0x80, // Some(-32759)
            0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [
            Some(&Value::Integer(-121)),
            Some(&Value::Integer(-120)),
            None,
        ];
        let expected = [
            0x12, // Some(Type::Int16(1))
            0x87, 0xff, // Some(-121)
            0x88, 0xff, // Some(-120)
            0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [
            Some(&Value::Integer(-120)),
            Some(&Value::Integer(-119)),
            None,
        ];
        let expected = [
            0x11, // Some(Type::Int8(1))
            0x88, // Some(-120)
            0x89, // Some(-119)
            0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [
            Some(&Value::Integer(-1)),
            Some(&Value::Integer(0)),
            Some(&Value::Integer(1)),
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

        let values = [Some(&Value::Integer(126)), Some(&Value::Integer(127)), None];
        let expected = [
            0x11, // Some(Type::Int8(1))
            0x7e, // Some(126)
            0x7f, // Some(127)
            0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [Some(&Value::Integer(127)), Some(&Value::Integer(128)), None];
        let expected = [
            0x12, // Some(Type::Int16(1))
            0x7f, 0x00, // Some(127)
            0x80, 0x00, // Some(128)
            0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let values = [
            Some(&Value::Integer(32766)),
            Some(&Value::Integer(32767)),
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
            Some(&Value::Integer(32767)),
            Some(&Value::Integer(32768)),
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
            Some(&Value::Integer(2147483646)),
            Some(&Value::Integer(2147483647)),
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
    fn test_write_genotype_field_values_with_integer_array_values() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            key: &Key,
            values: &[Option<&Value>],
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_genotype_field_values(buf, &key, &values)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let key = Key::Other(
            String::from("I32"),
            Number::Count(2),
            format::Type::Integer,
            String::default(),
        );

        let mut buf = Vec::new();

        let value_0 = Value::IntegerArray(vec![Some(-2147483641), Some(-2147483640)]);
        let value_1 = Value::IntegerArray(vec![Some(-2147483639), None]);
        let value_2 = Value::IntegerArray(vec![Some(-2147483638)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];
        buf.clear();
        assert!(matches!(
            write_genotype_field_values(&mut buf, &key, &values),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        let value_0 = Value::IntegerArray(vec![Some(-2147483640), Some(-2147483639)]);
        let value_1 = Value::IntegerArray(vec![Some(-2147483638), None]);
        let value_2 = Value::IntegerArray(vec![Some(-2147483637)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];
        #[rustfmt::skip]
        let expected = [
            0x23, // Some(Type::Int32(2))
            0x08, 0x00, 0x00, 0x80, 0x09, 0x00, 0x00, 0x80, // Some([Some(-2147483640), Some(-2147483639)])
            0x0a, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, // Some([Some(-2147483638), None])
            0x0b, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // Some([Some(-2147483637)])
            0x00, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let value_0 = Value::IntegerArray(vec![Some(-32761), Some(-32760)]);
        let value_1 = Value::IntegerArray(vec![Some(-32759), None]);
        let value_2 = Value::IntegerArray(vec![Some(-32758)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];
        #[rustfmt::skip]
        let expected = [
            0x23, // Some(Type::Int32(2))
            0x07, 0x80, 0xff, 0xff, 0x08, 0x80, 0xff, 0xff, // Some([Some(-32761), Some(-32760)])
            0x09, 0x80, 0xff, 0xff, 0x00, 0x00, 0x00, 0x80, // Some([Some(-32759), None])
            0x0a, 0x80, 0xff, 0xff, 0x01, 0x00, 0x00, 0x80, // Some([Some(-32761)])
            0x00, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let value_0 = Value::IntegerArray(vec![Some(-32760), Some(-32759)]);
        let value_1 = Value::IntegerArray(vec![Some(-32758), None]);
        let value_2 = Value::IntegerArray(vec![Some(-32757)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];
        let expected = [
            0x22, // Some(Type::Int16(2))
            0x08, 0x80, 0x09, 0x80, // Some([Some(-32760), Some(-32759)])
            0x0a, 0x80, 0x00, 0x80, // Some([Some(-32758), None])
            0x0b, 0x80, 0x01, 0x80, // Some([Some(-32757)])
            0x00, 0x80, 0x01, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let value_0 = Value::IntegerArray(vec![Some(-121), Some(-120)]);
        let value_1 = Value::IntegerArray(vec![Some(-119), None]);
        let value_2 = Value::IntegerArray(vec![Some(-118)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];
        let expected = [
            0x22, // Some(Type::Int16(2))
            0x87, 0xff, 0x88, 0xff, // Some([Some(-121), Some(-120)])
            0x89, 0xff, 0x00, 0x80, // Some([Some(-119), None])
            0x8a, 0xff, 0x01, 0x80, // Some([Some(-118)])
            0x00, 0x80, 0x01, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let value_0 = Value::IntegerArray(vec![Some(-120), Some(-119)]);
        let value_1 = Value::IntegerArray(vec![Some(-118), None]);
        let value_2 = Value::IntegerArray(vec![Some(-117)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];
        let expected = [
            0x21, // Some(Type::Int8(2))
            0x88, 0x89, // Some([Some(-120), Some(-119)])
            0x8a, 0x80, // Some([Some(-118), None])
            0x8b, 0x81, // Some([Some(-117)])
            0x80, 0x81, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let value_0 = Value::IntegerArray(vec![Some(-1), Some(0)]);
        let value_1 = Value::IntegerArray(vec![Some(1), None]);
        let value_2 = Value::IntegerArray(vec![Some(2)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];
        let expected = [
            0x21, // Some(Type::Int8(2))
            0xff, 0x00, // Some([Some(-1), Some(0)])
            0x01, 0x80, // Some([Some(1), None])
            0x02, 0x81, // Some([Some(2)])
            0x80, 0x81, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let value_0 = Value::IntegerArray(vec![Some(124), Some(125)]);
        let value_1 = Value::IntegerArray(vec![Some(126), None]);
        let value_2 = Value::IntegerArray(vec![Some(127)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];
        let expected = [
            0x21, // Some(Type::Int8(2))
            0x7c, 0x7d, // Some([Some(124), Some(125)])
            0x7e, 0x80, // Some([Some(126), None])
            0x7f, 0x81, // Some([Some(127)])
            0x80, 0x81, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let value_0 = Value::IntegerArray(vec![Some(125), Some(126)]);
        let value_1 = Value::IntegerArray(vec![Some(127), None]);
        let value_2 = Value::IntegerArray(vec![Some(128)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];
        let expected = [
            0x22, // Some(Type::Int16(2))
            0x7d, 0x00, 0x7e, 0x00, // Some([Some(125), Some(126)])
            0x7f, 0x00, 0x00, 0x80, // Some([Some(127), None])
            0x80, 0x00, 0x01, 0x80, // Some([Some(128)])
            0x00, 0x80, 0x01, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let value_0 = Value::IntegerArray(vec![Some(32764), Some(32765)]);
        let value_1 = Value::IntegerArray(vec![Some(32766), None]);
        let value_2 = Value::IntegerArray(vec![Some(32767)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];
        let expected = [
            0x22, // Some(Type::Int16(2))
            0xfc, 0x7f, 0xfd, 0x7f, // Some([Some(32764), Some(32765)])
            0xfe, 0x7f, 0x00, 0x80, // Some([Some(32766), None])
            0xff, 0x7f, 0x01, 0x80, // Some([Some(32767)])
            0x00, 0x80, 0x01, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let value_0 = Value::IntegerArray(vec![Some(32765), Some(32766)]);
        let value_1 = Value::IntegerArray(vec![Some(32767), None]);
        let value_2 = Value::IntegerArray(vec![Some(32768)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];
        #[rustfmt::skip]
        let expected = [
            0x23, // Some(Type::Int32(2))
            0xfd, 0x7f, 0x00, 0x00, 0xfe, 0x7f, 0x00, 0x00, // Some([Some(32765), Some(32766)])
            0xff, 0x7f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, // Some([Some(32767), None])
            0x00, 0x80, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, // Some([Some(32768)])
            0x00, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        let value_0 = Value::IntegerArray(vec![Some(2147483644), Some(2147483645)]);
        let value_1 = Value::IntegerArray(vec![Some(2147483646), None]);
        let value_2 = Value::IntegerArray(vec![Some(2147483647)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];
        #[rustfmt::skip]
        let expected = [
            0x23, // Some(Type::Int32(2))
            0xfc, 0xff, 0xff, 0x7f, 0xfd, 0xff, 0xff, 0x7f, // Some([Some(2147483645), Some(2147483645)])
            0xfe, 0xff, 0xff, 0x7f, 0x00, 0x00, 0x00, 0x80, // Some([Some(2147483646), None])
            0xff, 0xff, 0xff, 0x7f, 0x01, 0x00, 0x00, 0x80, // Some([Some(2147483647)])
            0x00, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // None
        ];
        t(&mut buf, &key, &values, &expected)?;

        Ok(())
    }

    #[test]
    fn test_write_genotype_field_values_with_float_values() -> io::Result<()> {
        let key = Key::Other(
            String::from("F32"),
            Number::Count(1),
            format::Type::Float,
            String::default(),
        );

        let values = [Some(&Value::Float(0.0)), Some(&Value::Float(1.0)), None];

        let mut buf = Vec::new();
        write_genotype_field_values(&mut buf, &key, &values)?;

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
    fn test_write_genotype_field_values_with_float_array_values() -> io::Result<()> {
        let key = Key::Other(
            String::from("F32"),
            Number::Count(2),
            format::Type::Float,
            String::default(),
        );

        let value_0 = Value::FloatArray(vec![Some(0.0), Some(1.0)]);
        let value_1 = Value::FloatArray(vec![Some(0.0), None]);
        let value_2 = Value::FloatArray(vec![Some(0.0)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];

        let mut buf = Vec::new();
        write_genotype_field_values(&mut buf, &key, &values)?;

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
    fn test_write_genotype_field_values_with_character_values() -> io::Result<()> {
        let key = Key::Other(
            String::from("CHAR"),
            Number::Count(1),
            format::Type::Character,
            String::default(),
        );

        let values = [
            Some(&Value::Character('n')),
            Some(&Value::Character('d')),
            None,
        ];

        let mut buf = Vec::new();
        write_genotype_field_values(&mut buf, &key, &values)?;

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
    fn test_write_genotype_field_values_with_character_array_values() -> io::Result<()> {
        let key = Key::Other(
            String::from("CHAR"),
            Number::Count(2),
            format::Type::Character,
            String::default(),
        );

        let value_0 = Value::CharacterArray(vec![Some('n'), Some('d')]);
        let value_1 = Value::CharacterArray(vec![Some('l'), None]);
        let value_2 = Value::CharacterArray(vec![Some('s')]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];

        let mut buf = Vec::new();
        write_genotype_field_values(&mut buf, &key, &values)?;

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
    fn test_write_genotype_field_values_with_string_values() -> io::Result<()> {
        let key = Key::Other(
            String::from("STRING"),
            Number::Count(1),
            format::Type::String,
            String::default(),
        );

        let value_0 = Value::String(String::from("n"));
        let value_1 = Value::String(String::from("ndl"));
        let value_2 = Value::String(String::from("ndls"));
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];

        let mut buf = Vec::new();
        write_genotype_field_values(&mut buf, &key, &values)?;

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
    fn test_write_genotype_field_values_with_string_array_values() -> io::Result<()> {
        let key = Key::Other(
            String::from("STRING"),
            Number::Count(2),
            format::Type::String,
            String::default(),
        );

        let value_0 = Value::StringArray(vec![Some(String::from("n")), Some(String::from("nd"))]);
        let value_1 = Value::StringArray(vec![Some(String::from("ndls")), None]);
        let value_2 = Value::StringArray(vec![Some(String::from("nd"))]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2), None];

        let mut buf = Vec::new();
        write_genotype_field_values(&mut buf, &key, &values)?;

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
}
