use std::{
    convert::TryFrom,
    io::{self, Read},
    str,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::{
    self as vcf,
    record::genotypes::{
        genotype::field::{Key, Value},
        Genotype,
    },
};

const NUL: u8 = 0x00;

use crate::{
    header::StringMap,
    reader::{string_map::read_string_map_index, value::read_type},
    record::value::{Float, Int16, Int32, Int8, Type},
};

pub fn read_genotypes<R>(
    reader: &mut R,
    string_map: &StringMap,
    sample_count: usize,
    format_count: usize,
) -> io::Result<Vec<Genotype>>
where
    R: Read,
{
    use vcf::record::genotypes::genotype::Field;

    let mut genotypes = vec![Vec::new(); sample_count];

    for _ in 0..format_count {
        let key = read_genotype_field_key(reader, string_map)?;

        let values = if key == Key::Genotype {
            read_genotype_genotype_field_values(reader, sample_count)?
        } else {
            read_genotype_field_values(reader, sample_count)?
        };

        for (fields, value) in genotypes.iter_mut().zip(values) {
            let field = Field::new(key.clone(), value);
            fields.push(field);
        }
    }

    genotypes
        .into_iter()
        .map(|fields| {
            Genotype::try_from(fields).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
        .collect()
}

fn read_genotype_field_key<R>(reader: &mut R, string_map: &StringMap) -> io::Result<Key>
where
    R: Read,
{
    read_string_map_index(reader)
        .and_then(|j| {
            string_map.get_index(j).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid string map index: {}", j),
                )
            })
        })
        .and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}

fn read_genotype_field_values<R>(
    reader: &mut R,
    sample_count: usize,
) -> io::Result<Vec<Option<Value>>>
where
    R: Read,
{
    match read_type(reader)? {
        Some(Type::Int8(len)) => match len {
            0 => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid number: {}", len),
            )),
            1 => read_genotype_field_int8_values(reader, sample_count),
            _ => read_genotype_field_int8_array_values(reader, sample_count, len),
        },
        Some(Type::Int16(len)) => match len {
            0 => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid number: {}", len),
            )),
            1 => read_genotype_field_int16_values(reader, sample_count),
            _ => read_genotype_field_int16_array_values(reader, sample_count, len),
        },
        Some(Type::Int32(len)) => match len {
            0 => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid number: {}", len),
            )),
            1 => read_genotype_field_int32_values(reader, sample_count),
            _ => read_genotype_field_int32_array_values(reader, sample_count, len),
        },
        Some(Type::Float(len)) => match len {
            0 => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid number: {}", len),
            )),
            1 => read_genotype_field_float_values(reader, sample_count),
            _ => read_genotype_field_float_array_values(reader, sample_count, len),
        },
        Some(Type::String(len)) => read_genotype_field_string_values(reader, sample_count, len),
        ty => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unhandled type: {:?}", ty),
        )),
    }
}

fn read_genotype_field_int8_values<R>(
    reader: &mut R,
    sample_count: usize,
) -> io::Result<Vec<Option<Value>>>
where
    R: Read,
{
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = reader.read_i8().map(Int8::from)?;

        match value {
            Int8::Value(n) => values.push(Some(Value::Integer(i32::from(n)))),
            Int8::Missing => values.push(None),
            _ => todo!("unhandled i8 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_genotype_field_int8_array_values<R>(
    reader: &mut R,
    sample_count: usize,
    len: usize,
) -> io::Result<Vec<Option<Value>>>
where
    R: Read,
{
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let mut buf = vec![0; len];
        reader.read_i8_into(&mut buf)?;

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
            values.push(Some(Value::IntegerArray(vs)));
        }
    }

    Ok(values)
}

fn read_genotype_field_int16_values<R>(
    reader: &mut R,
    sample_count: usize,
) -> io::Result<Vec<Option<Value>>>
where
    R: Read,
{
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = reader.read_i16::<LittleEndian>().map(Int16::from)?;

        match value {
            Int16::Value(n) => values.push(Some(Value::Integer(i32::from(n)))),
            Int16::Missing => values.push(None),
            _ => todo!("unhandled i16 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_genotype_field_int16_array_values<R>(
    reader: &mut R,
    sample_count: usize,
    len: usize,
) -> io::Result<Vec<Option<Value>>>
where
    R: Read,
{
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let mut buf = vec![0; len];
        reader.read_i16_into::<LittleEndian>(&mut buf)?;

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
            values.push(Some(Value::IntegerArray(vs)));
        }
    }

    Ok(values)
}

fn read_genotype_field_int32_values<R>(
    reader: &mut R,
    sample_count: usize,
) -> io::Result<Vec<Option<Value>>>
where
    R: Read,
{
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = reader.read_i32::<LittleEndian>().map(Int32::from)?;

        match value {
            Int32::Value(n) => values.push(Some(Value::Integer(n))),
            Int32::Missing => values.push(None),
            _ => todo!("unhandled i32 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_genotype_field_int32_array_values<R>(
    reader: &mut R,
    sample_count: usize,
    len: usize,
) -> io::Result<Vec<Option<Value>>>
where
    R: Read,
{
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let mut buf = vec![0; len];
        reader.read_i32_into::<LittleEndian>(&mut buf)?;

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
            values.push(Some(Value::IntegerArray(vs)));
        }
    }

    Ok(values)
}

fn read_genotype_field_float_values<R>(
    reader: &mut R,
    sample_count: usize,
) -> io::Result<Vec<Option<Value>>>
where
    R: Read,
{
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = reader.read_f32::<LittleEndian>().map(Float::from)?;

        match value {
            Float::Value(n) => values.push(Some(Value::Float(n))),
            Float::Missing => values.push(None),
            _ => todo!("unhandled f32 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_genotype_field_float_array_values<R>(
    reader: &mut R,
    sample_count: usize,
    len: usize,
) -> io::Result<Vec<Option<Value>>>
where
    R: Read,
{
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let mut buf = vec![0.0; len];
        reader.read_f32_into::<LittleEndian>(&mut buf)?;

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
            values.push(Some(Value::FloatArray(vs)));
        }
    }

    Ok(values)
}

fn read_genotype_field_string_values<R>(
    reader: &mut R,
    sample_count: usize,
    len: usize,
) -> io::Result<Vec<Option<Value>>>
where
    R: Read,
{
    let mut values = Vec::with_capacity(sample_count);
    let mut buf = vec![0; len];

    for _ in 0..sample_count {
        reader.read_exact(&mut buf)?;

        let data = match buf.iter().position(|&b| b == NUL) {
            Some(i) => &buf[..i],
            None => &buf[..],
        };

        let s = str::from_utf8(data).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let value = Value::String(s.into());

        values.push(Some(value));
    }

    Ok(values)
}

fn read_genotype_genotype_field_values<R>(
    reader: &mut R,
    sample_count: usize,
) -> io::Result<Vec<Option<Value>>>
where
    R: Read,
{
    let mut values = Vec::with_capacity(sample_count);

    match read_type(reader)? {
        Some(Type::Int8(len)) => match len {
            0 => values.push(None),
            1 => {
                for _ in 0..sample_count {
                    let value = reader
                        .read_i8()
                        .map(|v| parse_genotype_genotype_field_values(&[v]))
                        .map(Value::String)?;

                    values.push(Some(value));
                }
            }
            _ => {
                for _ in 0..sample_count {
                    let mut buf = vec![0; len];
                    reader.read_i8_into(&mut buf)?;
                    let value = Value::String(parse_genotype_genotype_field_values(&buf));
                    values.push(Some(value));
                }
            }
        },
        ty => todo!("unhandled type: {:?}", ty),
    }

    Ok(values)
}

fn parse_genotype_genotype_field_values(values: &[i8]) -> String {
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
            genotype.push_str(&format!("{}", j));
        }
    }

    genotype
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_genotype_field_values_with_int8_values() -> io::Result<()> {
        let data = [
            0x11, // Some(Type::Int8(1))
            0x05, // Some(5)
            0x08, // Some(8)
            0x80, // None
        ];
        let mut reader = &data[..];

        let actual = read_genotype_field_values(&mut reader, 3)?;
        let expected = vec![Some(Value::Integer(5)), Some(Value::Integer(8)), None];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_genotype_field_values_with_int8_array_values() -> io::Result<()> {
        let data = [
            0x21, // Some(Type::Int8(2))
            0x05, 0x08, // Some([Some(5), Some(8)])
            0x0d, 0x80, // Some([Some(13), None])
            0x15, 0x81, // Some([Some(21)])
            0x80, 0x81, // None
        ];
        let mut reader = &data[..];

        let actual = read_genotype_field_values(&mut reader, 4)?;
        let expected = vec![
            Some(Value::IntegerArray(vec![Some(5), Some(8)])),
            Some(Value::IntegerArray(vec![Some(13), None])),
            Some(Value::IntegerArray(vec![Some(21)])),
            None,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_genotype_field_values_with_int16_values() -> io::Result<()> {
        let data = [
            0x12, // Some(Type::Int16(1))
            0x05, 0x00, // Some(5)
            0x08, 0x00, // Some(8)
            0x00, 0x80, // None
        ];
        let mut reader = &data[..];

        let actual = read_genotype_field_values(&mut reader, 3)?;
        let expected = vec![Some(Value::Integer(5)), Some(Value::Integer(8)), None];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_genotype_field_values_with_int16_array_values() -> io::Result<()> {
        let data = [
            0x22, // Some(Type::Int16(2))
            0x05, 0x00, 0x08, 0x00, // Some([Some(5), Some(8)])
            0x0d, 0x00, 0x00, 0x80, // Some([Some(13), None])
            0x15, 0x00, 0x01, 0x80, // Some([Some(21)])
            0x00, 0x80, 0x01, 0x80, // None
        ];
        let mut reader = &data[..];

        let actual = read_genotype_field_values(&mut reader, 4)?;
        let expected = vec![
            Some(Value::IntegerArray(vec![Some(5), Some(8)])),
            Some(Value::IntegerArray(vec![Some(13), None])),
            Some(Value::IntegerArray(vec![Some(21)])),
            None,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_genotype_field_values_with_int32_values() -> io::Result<()> {
        let data = [
            0x13, // Some(Type::Int32(1))
            0x05, 0x00, 0x00, 0x00, // Some(5)
            0x08, 0x00, 0x00, 0x00, // Some(8)
            0x00, 0x00, 0x00, 0x80, // None
        ];
        let mut reader = &data[..];

        let actual = read_genotype_field_values(&mut reader, 3)?;
        let expected = vec![Some(Value::Integer(5)), Some(Value::Integer(8)), None];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_genotype_field_values_with_int32_array_values() -> io::Result<()> {
        let data = [
            0x23, // Some(Type::Int32(2))
            0x05, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x00, // Some([Some(5), Some(8)])
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, // Some([Some(13), None])
            0x15, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, // Some([Some(21)])
            0x00, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // None
        ];
        let mut reader = &data[..];

        let actual = read_genotype_field_values(&mut reader, 4)?;
        let expected = vec![
            Some(Value::IntegerArray(vec![Some(5), Some(8)])),
            Some(Value::IntegerArray(vec![Some(13), None])),
            Some(Value::IntegerArray(vec![Some(21)])),
            None,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_genotype_field_values_with_float_values() -> io::Result<()> {
        let data = [
            0x15, // Some(Type::Float(1))
            0x00, 0x00, 0x00, 0x00, // Some(0.0)
            0x00, 0x00, 0x80, 0x3f, // Some(1.0)
            0x01, 0x00, 0x80, 0x7f, // None
        ];
        let mut reader = &data[..];

        let actual = read_genotype_field_values(&mut reader, 3)?;
        let expected = vec![Some(Value::Float(0.0)), Some(Value::Float(1.0)), None];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_genotype_field_values_with_float_array_values() -> io::Result<()> {
        let data = [
            0x25, // Some(Type::Float(2))
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x3f, // Some([Some(0.0), Some(1.0)])
            0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x80, 0x7f, // Some([Some(0.0), None])
            0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x80, 0x7f, // Some([Some(0.0)])
            0x01, 0x00, 0x80, 0x7f, 0x02, 0x00, 0x80, 0x7f, // None
        ];
        let mut reader = &data[..];

        let actual = read_genotype_field_values(&mut reader, 4)?;
        let expected = vec![
            Some(Value::FloatArray(vec![Some(0.0), Some(1.0)])),
            Some(Value::FloatArray(vec![Some(0.0), None])),
            Some(Value::FloatArray(vec![Some(0.0)])),
            None,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_genotype_field_values_with_string_values() -> io::Result<()> {
        let data = [
            0x47, // Some(Type::String(4))
            b'n', 0x00, 0x00, 0x00, // "n"
            b'n', b'd', b'l', 0x00, // "ndl"
            b'n', b'd', b'l', b's', // "ndls"
        ];
        let mut reader = &data[..];

        let actual = read_genotype_field_values(&mut reader, 3)?;
        let expected = vec![
            Some(Value::String(String::from("n"))),
            Some(Value::String(String::from("ndl"))),
            Some(Value::String(String::from("ndls"))),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_parse_genotype_genotype_field_values() {
        // Examples from § 6.3.3 Type encoding (2021-05-13)

        assert_eq!(parse_genotype_genotype_field_values(&[0x02, 0x02]), "0/0");
        assert_eq!(parse_genotype_genotype_field_values(&[0x02, 0x04]), "0/1");
        assert_eq!(parse_genotype_genotype_field_values(&[0x04, 0x04]), "1/1");
        assert_eq!(parse_genotype_genotype_field_values(&[0x02, 0x05]), "0|1");
        assert_eq!(parse_genotype_genotype_field_values(&[0x00, 0x00]), "./.");
        assert_eq!(parse_genotype_genotype_field_values(&[0x02]), "0");
        assert_eq!(parse_genotype_genotype_field_values(&[0x04]), "1");
        assert_eq!(
            parse_genotype_genotype_field_values(&[0x02, 0x04, 0x06]),
            "0/1/2"
        );
        assert_eq!(
            parse_genotype_genotype_field_values(&[0x02, 0x04, 0x07]),
            "0/1|2"
        );
        assert_eq!(
            parse_genotype_genotype_field_values(&[0x02, i8::from(Int8::EndOfVector)]),
            "0"
        );
    }
}
