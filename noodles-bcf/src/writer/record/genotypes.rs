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
        write_genotype_key(writer, string_map, key)?;

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

        write_genotype_values(writer, key, &values)?;
    }

    Ok(())
}

pub fn write_genotype_key<W>(
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

pub fn write_genotype_values<W>(
    writer: &mut W,
    key: &vcf::record::genotype::field::Key,
    values: &[Option<&vcf::record::genotype::field::Value>],
) -> io::Result<()>
where
    W: Write,
{
    use vcf::{header::format, record::genotype};

    match key.ty() {
        format::Type::Integer => {
            write_type(writer, Some(Type::Int32(1)))?;

            for value in values {
                match value {
                    Some(genotype::field::Value::Integer(n)) => {
                        writer.write_i32::<LittleEndian>(*n)?;
                    }
                    Some(v) => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("expected integer value, got {:?}", v),
                        ))
                    }
                    None => writer.write_i32::<LittleEndian>(i32::from(Int32::Missing))?,
                }
            }
        }
        ty => todo!("unhandled genotype value: {:?}", ty),
    }

    Ok(())
}
