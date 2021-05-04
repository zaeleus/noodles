use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::ReadBytesExt;
use noodles_vcf::{self as vcf, record::Genotype};

use crate::{
    header::StringMap,
    reader::value::{read_type, read_value, Value},
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
    use crate::reader::value::Int8;

    match read_value(reader)? {
        Some(Value::Int8(Some(Int8::Value(i)))) => usize::try_from(i)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
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
            }),
        v => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected i8, got {:?}", v),
            ));
        }
    }
}

fn read_genotype_values<R>(
    reader: &mut R,
    sample_count: usize,
) -> io::Result<Vec<Option<vcf::record::genotype::field::Value>>>
where
    R: Read,
{
    use vcf::record::genotype;

    use crate::reader::value::{Int8, Type};

    let mut values = Vec::with_capacity(sample_count);

    match read_type(reader)? {
        Some(Type::Int8(len)) => match len {
            0 => values.push(None),
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
                                _ => todo!("unhanlded i8 array value: {:?}", value),
                            })
                            .collect(),
                    );
                    values.push(Some(value));
                }
            }
        },
        ty => todo!("unhandled type: {:?}", ty),
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

    use crate::reader::value::Type;

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
    use crate::reader::value::Int8;

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
    fn test_parse_genotype_genotype_values() {
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x02]), "0/0");
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x04]), "0/1");
        assert_eq!(parse_genotype_genotype_values(&[0x04, 0x04]), "1/1");
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x05]), "0|1");
        assert_eq!(parse_genotype_genotype_values(&[0x00, 0x00]), "./.");
        assert_eq!(parse_genotype_genotype_values(&[0x02]), "0");
        assert_eq!(parse_genotype_genotype_values(&[0x04]), "1");
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x04, 0x06]), "0/1/2");
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x04, 0x07]), "0/1|2");
        assert_eq!(parse_genotype_genotype_values(&[0x02, -127]), "0");
    }
}
