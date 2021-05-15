use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::{self as vcf, record::Genotype};

use crate::{
    header::StringMap,
    reader::{string_map::read_string_map_index, value::read_type},
    record::value::{Int16, Int32, Int8, Type},
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
    use vcf::record::genotype::{self, Field};

    let mut genotypes = vec![Vec::new(); sample_count];

    for _ in 0..format_count {
        let key = read_genotype_key(reader, string_map)?;

        let values = if key == genotype::field::Key::Genotype {
            read_genotype_genotype_values(reader, sample_count)?
        } else {
            read_genotype_values(reader, sample_count)?
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

fn read_genotype_key<R>(
    reader: &mut R,
    string_map: &StringMap,
) -> io::Result<vcf::record::genotype::field::Key>
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

fn read_genotype_values<R>(
    reader: &mut R,
    sample_count: usize,
) -> io::Result<Vec<Option<vcf::record::genotype::field::Value>>>
where
    R: Read,
{
    match read_type(reader)? {
        Some(Type::Int8(len)) => read_genotype_int8_values(reader, sample_count, len),
        Some(Type::Int16(len)) => read_genotype_int16_values(reader, sample_count, len),
        Some(Type::Int32(len)) => read_genotype_int32_values(reader, sample_count, len),
        ty => todo!("unhandled type: {:?}", ty),
    }
}

fn read_genotype_int8_values<R>(
    reader: &mut R,
    sample_count: usize,
    len: usize,
) -> io::Result<Vec<Option<vcf::record::genotype::field::Value>>>
where
    R: Read,
{
    use vcf::record::genotype;

    let mut values = Vec::with_capacity(sample_count);

    match len {
        0 => todo!("unhandled genotypes i8 len: 0"),
        1 => {
            for _ in 0..sample_count {
                let value = reader.read_i8().map(Int8::from)?;

                match value {
                    Int8::Value(n) => {
                        values.push(Some(genotype::field::Value::Integer(i32::from(n))))
                    }
                    Int8::Missing => values.push(None),
                    _ => todo!("unhandled i8 value: {:?}", value),
                }
            }
        }
        _ => {
            for _ in 0..sample_count {
                let mut buf = vec![0; len];
                reader.read_i8_into(&mut buf)?;

                let value = genotype::field::Value::IntegerArray(
                    buf.into_iter()
                        .map(Int8::from)
                        .map(|value| match value {
                            Int8::Value(n) => Some(i32::from(n)),
                            Int8::Missing => None,
                            _ => todo!("unhandled i8 array value: {:?}", value),
                        })
                        .collect(),
                );

                values.push(Some(value));
            }
        }
    }

    Ok(values)
}

fn read_genotype_int16_values<R>(
    reader: &mut R,
    sample_count: usize,
    len: usize,
) -> io::Result<Vec<Option<vcf::record::genotype::field::Value>>>
where
    R: Read,
{
    use vcf::record::genotype;

    let mut values = Vec::with_capacity(sample_count);

    match len {
        0 => todo!("unhandled genotypes i16 len: 0"),
        1 => {
            for _ in 0..sample_count {
                let value = reader.read_i16::<LittleEndian>().map(Int16::from)?;

                match value {
                    Int16::Value(n) => {
                        values.push(Some(genotype::field::Value::Integer(i32::from(n))))
                    }
                    Int16::Missing => values.push(None),
                    _ => todo!("unhandled i16 value: {:?}", value),
                }
            }
        }
        _ => {
            for _ in 0..sample_count {
                let mut buf = vec![0; len];
                reader.read_i16_into::<LittleEndian>(&mut buf)?;

                let value = genotype::field::Value::IntegerArray(
                    buf.into_iter()
                        .map(Int16::from)
                        .map(|value| match value {
                            Int16::Value(n) => Some(i32::from(n)),
                            Int16::Missing => None,
                            _ => todo!("unhandled i16 array value: {:?}", value),
                        })
                        .collect(),
                );

                values.push(Some(value));
            }
        }
    }

    Ok(values)
}

fn read_genotype_int32_values<R>(
    reader: &mut R,
    sample_count: usize,
    len: usize,
) -> io::Result<Vec<Option<vcf::record::genotype::field::Value>>>
where
    R: Read,
{
    use vcf::record::genotype;

    let mut values = Vec::with_capacity(sample_count);

    match len {
        0 => todo!("unhandled genotypes i32 len: 0"),
        1 => {
            for _ in 0..sample_count {
                let value = reader.read_i32::<LittleEndian>().map(Int32::from)?;

                match value {
                    Int32::Value(n) => values.push(Some(genotype::field::Value::Integer(n))),
                    Int32::Missing => values.push(None),
                    _ => todo!("unhandled i32 value: {:?}", value),
                }
            }
        }
        _ => {
            for _ in 0..sample_count {
                let mut buf = vec![0; len];
                reader.read_i32_into::<LittleEndian>(&mut buf)?;

                let value = genotype::field::Value::IntegerArray(
                    buf.into_iter()
                        .map(Int32::from)
                        .map(|value| match value {
                            Int32::Value(n) => Some(n),
                            Int32::Missing => None,
                            _ => todo!("unhandled i32 array value: {:?}", value),
                        })
                        .collect(),
                );

                values.push(Some(value));
            }
        }
    }

    Ok(values)
}

fn read_genotype_genotype_values<R>(
    reader: &mut R,
    sample_count: usize,
) -> io::Result<Vec<Option<vcf::record::genotype::field::Value>>>
where
    R: Read,
{
    use vcf::record::genotype;

    let mut values = Vec::with_capacity(sample_count);

    match read_type(reader)? {
        Some(Type::Int8(len)) => match len {
            0 => values.push(None),
            1 => {
                for _ in 0..sample_count {
                    let value = reader
                        .read_i8()
                        .map(|v| parse_genotype_genotype_values(&[v]))
                        .map(genotype::field::Value::String)?;

                    values.push(Some(value));
                }
            }
            _ => {
                for _ in 0..sample_count {
                    let mut buf = vec![0; len];
                    reader.read_i8_into(&mut buf)?;
                    let value =
                        genotype::field::Value::String(parse_genotype_genotype_values(&buf));
                    values.push(Some(value));
                }
            }
        },
        ty => todo!("unhandled type: {:?}", ty),
    }

    Ok(values)
}

fn parse_genotype_genotype_values(values: &[i8]) -> String {
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
    use vcf::record::genotype::field::Value;

    use super::*;

    const SAMPLE_COUNT: usize = 3;

    #[test]
    fn test_read_genotype_values_with_int8_values() -> io::Result<()> {
        let data = [
            0x11, // Some(Type::Int8(1))
            0x05, // Some(5)
            0x08, // Some(8)
            0x80, // None
        ];
        let mut reader = &data[..];
        let actual = read_genotype_values(&mut reader, SAMPLE_COUNT)?;
        let expected = vec![Some(Value::Integer(5)), Some(Value::Integer(8)), None];
        assert_eq!(actual, expected);

        let data = [
            0x21, // Some(Type::Int8(2))
            0x05, 0x08, // [Some(5), Some(8)]
            0x0d, 0x15, // [Some(13), Some(21)]
            0x80, 0x80, // [None, None]
        ];
        let mut reader = &data[..];
        let actual = read_genotype_values(&mut reader, SAMPLE_COUNT)?;
        let expected = vec![
            Some(Value::IntegerArray(vec![Some(5), Some(8)])),
            Some(Value::IntegerArray(vec![Some(13), Some(21)])),
            Some(Value::IntegerArray(vec![None, None])),
        ];
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_genotype_values_with_int16_values() -> io::Result<()> {
        let data = [
            0x12, // Some(Type::Int16(1))
            0x05, 0x00, // Some(5)
            0x08, 0x00, // Some(8)
            0x00, 0x80, // None
        ];
        let mut reader = &data[..];
        let actual = read_genotype_values(&mut reader, SAMPLE_COUNT)?;
        let expected = vec![Some(Value::Integer(5)), Some(Value::Integer(8)), None];
        assert_eq!(actual, expected);

        let data = [
            0x22, // Some(Type::Int16(2))
            0x05, 0x00, 0x08, 0x00, // [Some(5), Some(8)]
            0x0d, 0x00, 0x15, 0x00, // [Some(13), Some(21)]
            0x00, 0x80, 0x00, 0x80, // [None, None]
        ];
        let mut reader = &data[..];
        let actual = read_genotype_values(&mut reader, SAMPLE_COUNT)?;
        let expected = vec![
            Some(Value::IntegerArray(vec![Some(5), Some(8)])),
            Some(Value::IntegerArray(vec![Some(13), Some(21)])),
            Some(Value::IntegerArray(vec![None, None])),
        ];
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_genotype_values_with_int32_values() -> io::Result<()> {
        let data = [
            0x13, // Some(Type::Int32(1))
            0x05, 0x00, 0x00, 0x00, // Some(5)
            0x08, 0x00, 0x00, 0x00, // Some(8)
            0x00, 0x00, 0x00, 0x80, // None
        ];
        let mut reader = &data[..];
        let actual = read_genotype_values(&mut reader, SAMPLE_COUNT)?;
        let expected = vec![Some(Value::Integer(5)), Some(Value::Integer(8)), None];
        assert_eq!(actual, expected);

        let data = [
            0x23, // Some(Type::Int32(2))
            0x05, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x00, // [Some(5), Some(8)]
            0x0d, 0x00, 0x00, 0x00, 0x15, 0x00, 0x00, 0x00, // [Some(13), Some(21)]
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, // [None, None]
        ];
        let mut reader = &data[..];
        let actual = read_genotype_values(&mut reader, SAMPLE_COUNT)?;
        let expected = vec![
            Some(Value::IntegerArray(vec![Some(5), Some(8)])),
            Some(Value::IntegerArray(vec![Some(13), Some(21)])),
            Some(Value::IntegerArray(vec![None, None])),
        ];
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_parse_genotype_genotype_values() {
        use crate::record::value::Int8;

        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x02]), "0/0");
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x04]), "0/1");
        assert_eq!(parse_genotype_genotype_values(&[0x04, 0x04]), "1/1");
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x05]), "0|1");
        assert_eq!(parse_genotype_genotype_values(&[0x00, 0x00]), "./.");
        assert_eq!(parse_genotype_genotype_values(&[0x02]), "0");
        assert_eq!(parse_genotype_genotype_values(&[0x04]), "1");
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x04, 0x06]), "0/1/2");
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x04, 0x07]), "0/1|2");
        assert_eq!(
            parse_genotype_genotype_values(&[0x02, i8::from(Int8::EndOfVector)]),
            "0"
        );
    }
}
