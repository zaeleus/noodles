use std::{
    convert::TryFrom,
    io::{self, Read},
};

use noodles_vcf::{self as vcf, header::info::Type};

use crate::{
    header::StringMap,
    reader::{string_map::read_string_map_index, value::read_value},
    record::{
        value::{Float, Int16, Int32, Int8},
        Value,
    },
};

pub fn read_info<R>(
    reader: &mut R,
    infos: &vcf::header::Infos,
    string_map: &StringMap,
    len: usize,
) -> io::Result<vcf::record::Info>
where
    R: Read,
{
    let mut fields = Vec::with_capacity(len);

    for _ in 0..len {
        let key = read_info_field_key(reader, string_map)?;

        let info = infos.get(&key).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing header INFO record for {}", key),
            )
        })?;

        let value = read_info_field_value(reader, &info)?;

        let field = vcf::record::info::Field::new(key, value);
        fields.push(field);
    }

    vcf::record::Info::try_from(fields).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_info_field_key<R>(
    reader: &mut R,
    string_map: &StringMap,
) -> io::Result<vcf::record::info::field::Key>
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

fn read_info_field_value<R>(
    reader: &mut R,
    info: &vcf::header::Info,
) -> io::Result<vcf::record::info::field::Value>
where
    R: Read,
{
    match info.ty() {
        Type::Integer => read_info_field_integer_value(reader),
        Type::Flag => read_info_field_flag_value(reader),
        Type::Float => read_info_field_float_value(reader),
        Type::Character => read_info_field_character_value(reader),
        Type::String => read_info_field_string_value(reader),
    }
}

fn read_info_field_integer_value<R>(reader: &mut R) -> io::Result<vcf::record::info::field::Value>
where
    R: Read,
{
    let value = match read_value(reader)? {
        Some(Value::Int8(Some(Int8::Value(n)))) => {
            vcf::record::info::field::Value::Integer(i32::from(n))
        }
        Some(Value::Int8Array(values)) => vcf::record::info::field::Value::IntegerArray(
            values.into_iter().map(i32::from).collect(),
        ),
        Some(Value::Int16(Some(Int16::Value(n)))) => {
            vcf::record::info::field::Value::Integer(i32::from(n))
        }
        Some(Value::Int16Array(values)) => vcf::record::info::field::Value::IntegerArray(
            values.into_iter().map(i32::from).collect(),
        ),
        Some(Value::Int32(Some(Int32::Value(n)))) => vcf::record::info::field::Value::Integer(n),
        Some(Value::Int32Array(values)) => vcf::record::info::field::Value::IntegerArray(values),
        v => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("type mismatch: expected {}, got {:?}", Type::Integer, v),
            ));
        }
    };

    Ok(value)
}

fn read_info_field_flag_value<R>(reader: &mut R) -> io::Result<vcf::record::info::field::Value>
where
    R: Read,
{
    match read_value(reader)? {
        Some(Value::Int8(Some(Int8::Value(1)))) | None => Ok(vcf::record::info::field::Value::Flag),
        v => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("type mismatch: expected {}, got {:?}", Type::Flag, v),
        )),
    }
}

fn read_info_field_float_value<R>(reader: &mut R) -> io::Result<vcf::record::info::field::Value>
where
    R: Read,
{
    match read_value(reader)? {
        Some(Value::Float(Some(Float::Value(n)))) => Ok(vcf::record::info::field::Value::Float(n)),
        Some(Value::FloatArray(values)) => Ok(vcf::record::info::field::Value::FloatArray(values)),
        v => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("type mismatch: expected {}, got {:?}", Type::Float, v),
        )),
    }
}

fn read_info_field_character_value<R>(reader: &mut R) -> io::Result<vcf::record::info::field::Value>
where
    R: Read,
{
    match read_value(reader)? {
        Some(Value::String(Some(s))) => s
            .chars()
            .next()
            .map(vcf::record::info::field::Value::Character)
            .map(Ok)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "INFO character value missing")
            })?,
        v => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("type mismatch: expected {}, got {:?}", Type::Character, v),
        )),
    }
}

fn read_info_field_string_value<R>(reader: &mut R) -> io::Result<vcf::record::info::field::Value>
where
    R: Read,
{
    match read_value(reader)? {
        Some(Value::String(Some(s))) => Ok(vcf::record::info::field::Value::String(s)),
        v => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("type mismatch: expected {}, got {:?}", Type::String, v),
        )),
    }
}

#[cfg(test)]
mod tests {
    use vcf::{header::Number, record::info::field::Key};

    use super::*;

    #[test]
    fn test_read_info_field_value_with_integer_value() -> io::Result<()> {
        fn t(data: &[u8], info: &vcf::header::Info, expected_value: i32) -> io::Result<()> {
            let mut reader = data;

            let actual = read_info_field_value(&mut reader, &info)?;
            let expected = vcf::record::info::field::Value::Integer(expected_value);

            assert_eq!(actual, expected);

            Ok(())
        }

        let info = vcf::header::Info::from(Key::Other(
            String::from("I32"),
            Number::Count(1),
            Type::Integer,
            String::default(),
        ));

        // Some(Value::Int8(Some(Int8::Value(8))))
        t(&[0x11, 0x08], &info, 8)?;
        // Some(Value::Int16(Some(Int16::Value(13))))
        t(&[0x12, 0x0d, 0x00], &info, 13)?;
        // Some(Value::Int32(Some(Int32::Value(21))))
        t(&[0x13, 0x15, 0x00, 0x00, 0x00], &info, 21)?;

        Ok(())
    }

    #[test]
    fn test_read_info_field_value_with_flag_value() -> io::Result<()> {
        fn t(data: &[u8], info: &vcf::header::Info) -> io::Result<()> {
            let mut reader = data;

            let actual = read_info_field_value(&mut reader, &info)?;
            let expected = vcf::record::info::field::Value::Flag;

            assert_eq!(actual, expected);

            Ok(())
        }

        let info = vcf::header::Info::from(Key::Other(
            String::from("BOOL"),
            Number::Count(1),
            Type::Flag,
            String::default(),
        ));

        // None
        t(&[0x00], &info)?;
        // Some(Value::Int8(Some(Int8::Value(1))))
        t(&[0x11, 0x01], &info)?;

        Ok(())
    }

    #[test]
    fn test_read_info_field_value_with_float_value() -> io::Result<()> {
        // Some(Value::Float(Some(Float::Value(0.0))))
        let data = [0x15, 0x00, 0x00, 0x00, 0x00];
        let mut reader = &data[..];

        let info = vcf::header::Info::from(Key::Other(
            String::from("F32"),
            Number::Count(1),
            Type::Float,
            String::default(),
        ));

        let actual = read_info_field_value(&mut reader, &info)?;
        let expected = vcf::record::info::field::Value::Float(0.0);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_info_field_value_with_character_value() -> io::Result<()> {
        // Some(Value::String(Some(String::from("n"))))
        let data = [0x17, 0x6e];
        let mut reader = &data[..];

        let info = vcf::header::Info::from(Key::Other(
            String::from("CHAR"),
            Number::Count(1),
            Type::Character,
            String::default(),
        ));

        let actual = read_info_field_value(&mut reader, &info)?;
        let expected = vcf::record::info::field::Value::Character('n');

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_info_field_value_with_string_value() -> io::Result<()> {
        // Some(Value::String(Some(String::from("ndls"))))
        let data = [0x47, 0x6e, 0x64, 0x6c, 0x73];
        let mut reader = &data[..];

        let info = vcf::header::Info::from(Key::Other(
            String::from("STRING"),
            Number::Count(1),
            Type::String,
            String::default(),
        ));

        let actual = read_info_field_value(&mut reader, &info)?;
        let expected = vcf::record::info::field::Value::String(String::from("ndls"));

        assert_eq!(actual, expected);

        Ok(())
    }
}
