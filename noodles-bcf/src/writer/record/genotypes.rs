use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_vcf as vcf;

use crate::{
    header::StringMap,
    record::value::{Int32, Type},
    writer::{string_map::write_string_map_index, value::write_type},
};

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
    key: &vcf::record::genotype::field::Key,
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
    key: &vcf::record::genotype::field::Key,
    values: &[Option<&vcf::record::genotype::field::Value>],
) -> io::Result<()>
where
    W: Write,
{
    use vcf::header::{format, Number};

    match key.ty() {
        format::Type::Integer => match key.number() {
            Number::Count(1) => write_genotype_field_integer_values(writer, values),
            n => todo!("unhandled genotype integer number: {:?}", n),
        },
        ty => todo!("unhandled genotype value: {:?}", ty),
    }
}

fn write_genotype_field_integer_values<W>(
    writer: &mut W,
    values: &[Option<&vcf::record::genotype::field::Value>],
) -> io::Result<()>
where
    W: Write,
{
    use vcf::record::genotype;

    write_type(writer, Some(Type::Int32(1)))?;

    for value in values {
        match value {
            Some(genotype::field::Value::Integer(n)) => {
                writer.write_i32::<LittleEndian>(*n)?;
            }
            Some(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected Integer, got {:?}", v),
                ))
            }
            None => writer.write_i32::<LittleEndian>(i32::from(Int32::Missing))?,
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{
        header::{format, Number},
        record::genotype::{self, field::Key},
    };

    use super::*;

    #[test]
    fn test_write_genotype_field_values_with_integer_values() -> io::Result<()> {
        let key = Key::Other(
            String::from("I32"),
            Number::Count(1),
            format::Type::Integer,
            String::default(),
        );

        let values = [
            Some(&genotype::field::Value::Integer(5)),
            Some(&genotype::field::Value::Integer(8)),
            None,
        ];

        let mut buf = Vec::new();
        write_genotype_field_values(&mut buf, &key, &values)?;

        let expected = [
            0x13, // Some(Type::Int32(1))
            0x05, 0x00, 0x00, 0x00, // Some(5)
            0x08, 0x00, 0x00, 0x00, // Some(8)
            0x00, 0x00, 0x00, 0x80, // None
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
