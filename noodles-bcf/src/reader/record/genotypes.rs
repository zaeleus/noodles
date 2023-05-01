use std::{
    io::{self, Read},
    str,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::{
    self as vcf,
    record::{
        genotypes::{
            keys::{key, Key},
            sample::Value,
            Keys,
        },
        Genotypes,
    },
};

const NUL: u8 = 0x00;

use crate::{
    header::string_maps::StringStringMap,
    lazy::record::value::{Float, Int16, Int32, Int8, Type},
    reader::{string_map::read_string_map_index, value::read_type},
};

pub fn read_genotypes<R>(
    reader: &mut R,
    formats: &vcf::header::Formats,
    string_map: &StringStringMap,
    sample_count: usize,
    format_count: usize,
) -> io::Result<Genotypes>
where
    R: Read,
{
    let mut keys = Vec::with_capacity(format_count);
    let mut values = vec![Vec::new(); sample_count];

    for _ in 0..format_count {
        let key = read_genotype_field_key(reader, formats, string_map)?;

        let vs = if key == key::GENOTYPE {
            read_genotype_genotype_field_values(reader, sample_count)?
        } else {
            read_genotype_field_values(reader, sample_count)?
        };

        keys.push(key);

        for (sample, value) in values.iter_mut().zip(vs) {
            sample.push(value);
        }
    }

    let keys = Keys::try_from(keys).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(Genotypes::new(keys, values))
}

fn read_genotype_field_key<R>(
    reader: &mut R,
    formats: &vcf::header::Formats,
    string_map: &StringStringMap,
) -> io::Result<Key>
where
    R: Read,
{
    read_string_map_index(reader)
        .and_then(|j| {
            string_map.get_index(j).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid string map index: {j}"),
                )
            })
        })
        .and_then(|raw_key| {
            formats
                .keys()
                .find(|k| k.as_ref() == raw_key)
                .cloned()
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("missing header FORMAT record for {raw_key}"),
                    )
                })
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
                format!("invalid number: {len}"),
            )),
            1 => read_genotype_field_int8_values(reader, sample_count),
            _ => read_genotype_field_int8_array_values(reader, sample_count, len),
        },
        Some(Type::Int16(len)) => match len {
            0 => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid number: {len}"),
            )),
            1 => read_genotype_field_int16_values(reader, sample_count),
            _ => read_genotype_field_int16_array_values(reader, sample_count, len),
        },
        Some(Type::Int32(len)) => match len {
            0 => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid number: {len}"),
            )),
            1 => read_genotype_field_int32_values(reader, sample_count),
            _ => read_genotype_field_int32_array_values(reader, sample_count, len),
        },
        Some(Type::Float(len)) => match len {
            0 => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid number: {len}"),
            )),
            1 => read_genotype_field_float_values(reader, sample_count),
            _ => read_genotype_field_float_array_values(reader, sample_count, len),
        },
        Some(Type::String(len)) => read_genotype_field_string_values(reader, sample_count, len),
        ty => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unhandled type: {ty:?}"),
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
            Int8::Value(n) => values.push(Some(Value::from(i32::from(n)))),
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
            values.push(Some(Value::from(vs)));
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
            Int16::Value(n) => values.push(Some(Value::from(i32::from(n)))),
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
            values.push(Some(Value::from(vs)));
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
            Int32::Value(n) => values.push(Some(Value::from(n))),
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
            values.push(Some(Value::from(vs)));
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
            Float::Value(n) => values.push(Some(Value::from(n))),
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
            values.push(Some(Value::from(vs)));
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
        let value = Value::from(s);

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
                        .map(Value::from)?;

                    values.push(Some(value));
                }
            }
            _ => {
                let mut buf = vec![0; len];

                for _ in 0..sample_count {
                    reader.read_i8_into(&mut buf)?;
                    let value = Value::from(parse_genotype_genotype_field_values(&buf));
                    values.push(Some(value));
                }
            }
        },
        ty => todo!("unhandled type: {:?}", ty),
    }

    Ok(values)
}

fn parse_genotype_genotype_field_values(values: &[i8]) -> String {
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
        let expected = vec![Some(Value::from(5)), Some(Value::from(8)), None];

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
            Some(Value::from(vec![Some(5), Some(8)])),
            Some(Value::from(vec![Some(13), None])),
            Some(Value::from(vec![Some(21)])),
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
        let expected = vec![Some(Value::from(5)), Some(Value::from(8)), None];

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
            Some(Value::from(vec![Some(5), Some(8)])),
            Some(Value::from(vec![Some(13), None])),
            Some(Value::from(vec![Some(21)])),
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
        let expected = vec![Some(Value::from(5)), Some(Value::from(8)), None];

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
            Some(Value::from(vec![Some(5), Some(8)])),
            Some(Value::from(vec![Some(13), None])),
            Some(Value::from(vec![Some(21)])),
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
        let expected = vec![Some(Value::from(0.0)), Some(Value::from(1.0)), None];

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
            Some(Value::from(vec![Some(0.0), Some(1.0)])),
            Some(Value::from(vec![Some(0.0), None])),
            Some(Value::from(vec![Some(0.0)])),
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
            Some(Value::from("n")),
            Some(Value::from("ndl")),
            Some(Value::from("ndls")),
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
