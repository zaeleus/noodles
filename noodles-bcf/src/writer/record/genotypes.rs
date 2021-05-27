use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_vcf as vcf;

use crate::{
    header::StringMap,
    record::value::{Float, Int32, Type},
    writer::{string_map::write_string_map_index, value::write_type},
};

const NUL: u8 = 0x00;
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
            _ => write_genotype_field_integer_array_values(writer, values),
        },
        format::Type::Float => match key.number() {
            Number::Count(1) => write_genotype_field_float_values(writer, values),
            _ => write_genotype_field_float_array_values(writer, values),
        },
        format::Type::String => match key.number() {
            Number::Count(1) => write_genotype_field_string_values(writer, values),
            _ => write_genotype_field_string_array_values(writer, values),
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

fn write_genotype_field_integer_array_values<W>(
    writer: &mut W,
    values: &[Option<&vcf::record::genotype::field::Value>],
) -> io::Result<()>
where
    W: Write,
{
    use vcf::record::genotype;

    let max_len = values
        .iter()
        .flat_map(|value| match value {
            Some(genotype::field::Value::IntegerArray(vs)) => Some(vs.len()),
            _ => None,
        })
        .max()
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "missing IntegerArray values")
        })?;

    write_type(writer, Some(Type::Int32(max_len)))?;

    for value in values {
        match value {
            Some(genotype::field::Value::IntegerArray(vs)) => {
                for v in vs {
                    let raw_value = v.unwrap_or(i32::from(Int32::Missing));
                    writer.write_i32::<LittleEndian>(raw_value)?;
                }

                if vs.len() < max_len {
                    for _ in 0..(max_len - vs.len()) {
                        writer.write_i32::<LittleEndian>(i32::from(Int32::EndOfVector))?;
                    }
                }
            }
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected IntegerArray, got {:?}", v),
                ))
            }
        }
    }

    Ok(())
}

fn write_genotype_field_float_values<W>(
    writer: &mut W,
    values: &[Option<&vcf::record::genotype::field::Value>],
) -> io::Result<()>
where
    W: Write,
{
    use vcf::record::genotype;

    write_type(writer, Some(Type::Float(1)))?;

    for value in values {
        match value {
            Some(genotype::field::Value::Float(n)) => {
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
    values: &[Option<&vcf::record::genotype::field::Value>],
) -> io::Result<()>
where
    W: Write,
{
    use vcf::record::genotype;

    let max_len = values
        .iter()
        .flat_map(|value| match value {
            Some(genotype::field::Value::FloatArray(vs)) => Some(vs.len()),
            _ => None,
        })
        .max()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing FloatArray values"))?;

    write_type(writer, Some(Type::Float(max_len)))?;

    for value in values {
        match value {
            Some(genotype::field::Value::FloatArray(vs)) => {
                for v in vs {
                    let raw_value = v.unwrap_or(f32::from(Float::Missing));
                    writer.write_f32::<LittleEndian>(raw_value)?;
                }

                if vs.len() < max_len {
                    for _ in 0..(max_len - vs.len()) {
                        writer.write_f32::<LittleEndian>(f32::from(Float::EndOfVector))?;
                    }
                }
            }
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("type mismatch: expected FloatArray, got {:?}", v),
                ))
            }
        }
    }

    Ok(())
}

fn write_genotype_field_string_values<W>(
    writer: &mut W,
    values: &[Option<&vcf::record::genotype::field::Value>],
) -> io::Result<()>
where
    W: Write,
{
    use vcf::record::genotype;

    let max_len = values
        .iter()
        .flat_map(|value| match value {
            Some(genotype::field::Value::String(s)) => Some(s.len()),
            _ => None,
        })
        .max()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing String values"))?;

    let mut buf = Vec::with_capacity(values.len() * max_len);

    for value in values {
        match value {
            Some(genotype::field::Value::String(s)) => {
                buf.extend(s.bytes());

                if s.len() < max_len {
                    buf.resize(buf.len() + (max_len - s.len()), NUL);
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
    values: &[Option<&vcf::record::genotype::field::Value>],
) -> io::Result<()>
where
    W: Write,
{
    use vcf::record::genotype;

    const DELIMITER: char = ',';

    let mut serialized_values = Vec::with_capacity(values.len());

    for value in values {
        match value {
            Some(genotype::field::Value::StringArray(vs)) => {
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

    #[test]
    fn test_write_genotype_field_values_with_integer_array_values() -> io::Result<()> {
        let key = Key::Other(
            String::from("I32"),
            Number::Count(2),
            format::Type::Integer,
            String::default(),
        );

        let value_0 = genotype::field::Value::IntegerArray(vec![Some(5), Some(8)]);
        let value_1 = genotype::field::Value::IntegerArray(vec![Some(13), None]);
        let value_2 = genotype::field::Value::IntegerArray(vec![Some(21)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2)];

        let mut buf = Vec::new();
        write_genotype_field_values(&mut buf, &key, &values)?;

        let expected = [
            0x23, // Some(Type::Int32(2))
            0x05, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x00, // [Some(5), Some(8)]
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, // [Some(13), None]
            0x15, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, // [Some(21)]
        ];

        assert_eq!(buf, expected);

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

        let values = [
            Some(&genotype::field::Value::Float(0.0)),
            Some(&genotype::field::Value::Float(1.0)),
            None,
        ];

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

        let value_0 = genotype::field::Value::FloatArray(vec![Some(0.0), Some(1.0)]);
        let value_1 = genotype::field::Value::FloatArray(vec![Some(0.0), None]);
        let value_2 = genotype::field::Value::FloatArray(vec![Some(0.0)]);
        let values = [Some(&value_0), Some(&value_1), Some(&value_2)];

        let mut buf = Vec::new();
        write_genotype_field_values(&mut buf, &key, &values)?;

        let expected = [
            0x25, // Some(Type::Float(2))
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x3f, // [Some(0.0), Some(1.0)]
            0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x80, 0x7f, // [Some(0.0), None]
            0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x80, 0x7f, // [Some(0.0)]
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

        let value_0 = genotype::field::Value::String(String::from("n"));
        let value_1 = genotype::field::Value::String(String::from("ndl"));
        let value_2 = genotype::field::Value::String(String::from("ndls"));
        let values = [Some(&value_0), Some(&value_1), Some(&value_2)];

        let mut buf = Vec::new();
        write_genotype_field_values(&mut buf, &key, &values)?;

        let expected = [
            0x47, // Some(Type::String(4))
            b'n', 0x00, 0x00, 0x00, // "n"
            b'n', b'd', b'l', 0x00, // "ndl"
            b'n', b'd', b'l', b's', // "ndls"
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

        let value_0 = genotype::field::Value::StringArray(vec![
            Some(String::from("n")),
            Some(String::from("nd")),
        ]);
        let value_1 = genotype::field::Value::StringArray(vec![Some(String::from("ndls")), None]);
        let value_2 = genotype::field::Value::StringArray(vec![Some(String::from("nd"))]);
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
